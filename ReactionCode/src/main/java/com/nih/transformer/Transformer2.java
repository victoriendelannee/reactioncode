package com.nih.transformer;

import static org.openscience.cdk.silent.SilentChemObjectBuilder.getInstance;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.silent.Bond;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.CircularFingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.tools.ExtAtomContainerManipulator;

import com.nih.fragments.FragmentUtils;
import com.nih.tools.tools;


/**
 * DEPRECATED
 * OLD Transformer version
 *
 */
public class Transformer2 {

	int stoichiometry = 5;
	boolean checkValence = false;
	
    /**
     * Apply a transform and get an unique product
     * @param reagents
     * @param tranform
     * @return
     * @throws CloneNotSupportedException
     * @throws CDKException
     */
    public IReaction applyTranform(IAtomContainerSet reagents, IReaction tranform) throws CloneNotSupportedException, CDKException {    	
    	IReaction reaction = SilentChemObjectBuilder.getInstance().newInstance(IReaction.class);
    	IAtomContainerSet smirksReactants = tranform.getReactants();
    	IAtomContainerSet smirksProducts = tranform.getProducts();
    	IAtomContainer aggregateProducts = SilentChemObjectBuilder.getInstance().newAtomContainer();
    	IAtomContainerSet reagentsCopy = SilentChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
    	
    	//prepare products container to apply the reaction and  and index bonds
    	int cpt = 1;
    	for (IAtomContainer reactant : reagents.atomContainers()) {
    		reagentsCopy.addAtomContainer(reactant.clone());
    		aggregateProducts.add(reactant);
    		for (int i = 0; i < reactant.getAtomCount(); i++) {
    			IAtom atom = reactant.getAtom(i);
    			atom.setID(cpt+"");
    			atom.setProperty(CDKConstants.ATOM_ATOM_MAPPING, cpt);
    			cpt++;
    		}
    	}

    	Map<IAtomContainer, int[]> mappings = new HashMap<IAtomContainer, int[]>();
    	//detect if all reactant patterns match
    	for (IAtomContainer smirksReactant : smirksReactants.atomContainers()) {
    		Pattern ptrn = Pattern.findSubstructure(smirksReactant);
    		int[] mapping = ptrn.match(aggregateProducts);
    		if (mapping.length > 0) {
    			mappings.put(smirksReactant, mapping);
    		}
    	}
    	
    	
    	//apply the reaction
    	if (!mappings.isEmpty()) {
    		//index atom product
    		Map<Integer,IAtom> indexProduct = new HashMap<Integer,IAtom>();
    		IAtomContainer aggregateSmirksProducts = SilentChemObjectBuilder.getInstance().newAtomContainer();
    		for (IAtomContainer smirksProduct : smirksProducts.atomContainers()) {
    			aggregateSmirksProducts.add(smirksProduct);
    			for (IAtom atom : smirksProduct.atoms()) {
    				indexProduct.put(atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING), atom);
    			}
    		}
    		
    		//apply changes
    		List<IBond> processed = new ArrayList<IBond>();
    		Map<IBond,IBond> changes = new HashMap<IBond,IBond>();
    		Map<Integer,IAtom> aam = new HashMap<Integer,IAtom>();
    		for (Entry<IAtomContainer, int[]> e : mappings.entrySet()) {
    			IAtomContainer smirksReactant = e.getKey();
    			int[] mapping = e.getValue();
    			for (int i = 0; i < mapping.length; i++) {
    				IAtom smirksReactantAtom = smirksReactant.getAtom(i);
    				IAtom smirksProductAtom = indexProduct.get(smirksReactantAtom.getProperty(CDKConstants.ATOM_ATOM_MAPPING));
    				aam.put(smirksReactantAtom.getProperty(CDKConstants.ATOM_ATOM_MAPPING), aggregateProducts.getAtom(mapping[i]));
    				    				
    				List<IBond> conSmirksReactants = smirksReactant.getConnectedBondsList(smirksReactantAtom);
    				List<IBond> conSmirksProducts = aggregateSmirksProducts.getConnectedBondsList(smirksProductAtom);
    				//remove process bond in reactants and products
    				for (IBond bond : new ArrayList<IBond>(conSmirksReactants)) {
    					if (processed.contains(bond)) {
    						conSmirksReactants.remove(bond);
    					}
    				}
    				for (IBond bond : new ArrayList<IBond>(conSmirksProducts)) {
    					if (processed.contains(bond)) {
    						conSmirksProducts.remove(bond);
    					}
    				}
    				
    				changes.putAll(findChanges(conSmirksReactants, conSmirksProducts));
    				
    				
    				processed.addAll(conSmirksReactants);
    				processed.addAll(conSmirksProducts);
    				
    			}
    			
    		}
			//apply bond and atom changes
    		applyChanges(aggregateProducts, changes,  aam, indexProduct);
    		//reset flags and get product separated
    		FragmentUtils.resetVisitedFlags(aggregateProducts);
    		IAtomContainerSet products = FragmentUtils.makeAtomContainerSet(aggregateProducts);
    		reaction.setReactants(reagentsCopy);
        	reaction.setProducts(products);

    	}
    	else {
    		System.err.println("Could not match the reaction");
    	}
		return reaction;
    	
    	//TODO check SMIRKS is correct (contains aam for all atoms)
    }

    
    /**
     * Apply a transform and get an unique product
     * @param reactantsSmiles
     * @param tranform
     * @return
     * @throws CloneNotSupportedException
     * @throws CDKException
     */
    public IReaction applyTranform(String reactantsSmiles, IReaction tranform) throws CloneNotSupportedException, CDKException {
    	IReaction reaction = SilentChemObjectBuilder.getInstance().newInstance(IReaction.class);
    	SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
    	String[] reactantsSmi = reactantsSmiles.split("\\.");
    	
    	IAtomContainerSet smirksReactants = tranform.getReactants();
    	IAtomContainerSet smirksProducts = tranform.getProducts();
    	IAtomContainer aggregateProducts = SilentChemObjectBuilder.getInstance().newAtomContainer();
    	IAtomContainerSet reagents = SilentChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
    	
    	//prepare products container to apply the reaction and  and index bonds
    	sp.kekulise(false);
    	int cpt = 1;
    	for (String reactant : reactantsSmi) {
    		IAtomContainer rea = sp.parseSmiles(reactant);
    		for (int i = 0; i < rea.getAtomCount(); i++) {
    			IAtom atom = rea.getAtom(i);
    			atom.setID(cpt+"");
    			atom.setProperty(CDKConstants.ATOM_ATOM_MAPPING, cpt);
    			cpt++;
    		}
    		aggregateProducts.add(rea);
    		reagents.addAtomContainer(rea.clone());
    	}

    	Map<IAtomContainer, int[]> mappings = new HashMap<IAtomContainer, int[]>();
    	//detect if all reactant patterns match
    	for (IAtomContainer smirksReactant : smirksReactants.atomContainers()) {
    		Pattern ptrn = Pattern.findSubstructure(smirksReactant);
    		int[] mapping = ptrn.match(aggregateProducts);
    		if (mapping.length > 0) {
    			mappings.put(smirksReactant, mapping);
    		}
    	}
    	
    	
    	//apply the reaction
    	if (!mappings.isEmpty()) {
    		//index atom product
    		Map<Integer,IAtom> indexProduct = new HashMap<Integer,IAtom>();
    		IAtomContainer aggregateSmirksProducts = SilentChemObjectBuilder.getInstance().newAtomContainer();
    		for (IAtomContainer smirksProduct : smirksProducts.atomContainers()) {
    			aggregateSmirksProducts.add(smirksProduct);
    			for (IAtom atom : smirksProduct.atoms()) {
    				indexProduct.put(atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING), atom);
    			}
    		}
    		
    		//apply changes
    		List<IBond> processed = new ArrayList<IBond>();
    		Map<IBond,IBond> changes = new HashMap<IBond,IBond>();
    		Map<Integer,IAtom> aam = new HashMap<Integer,IAtom>();
    		for (Entry<IAtomContainer, int[]> e : mappings.entrySet()) {
    			IAtomContainer smirksReactant = e.getKey();
    			int[] mapping = e.getValue();
    			for (int i = 0; i < mapping.length; i++) {
    				IAtom smirksReactantAtom = smirksReactant.getAtom(i);
    				IAtom smirksProductAtom = indexProduct.get(smirksReactantAtom.getProperty(CDKConstants.ATOM_ATOM_MAPPING));
    				aam.put(smirksReactantAtom.getProperty(CDKConstants.ATOM_ATOM_MAPPING), aggregateProducts.getAtom(mapping[i]));

    				List<IBond> conSmirksReactants = smirksReactant.getConnectedBondsList(smirksReactantAtom);
    				List<IBond> conSmirksProducts = aggregateSmirksProducts.getConnectedBondsList(smirksProductAtom);
    				//remove process bond in reactants and products
    				for (IBond bond : new ArrayList<IBond>(conSmirksReactants)) {
    					if (processed.contains(bond)) {
    						conSmirksReactants.remove(bond);
    					}
    				}
    				for (IBond bond : new ArrayList<IBond>(conSmirksProducts)) {
    					if (processed.contains(bond)) {
    						conSmirksProducts.remove(bond);
    					}
    				}
    				
    				changes.putAll(findChanges(conSmirksReactants, conSmirksProducts));
    				
    				
    				processed.addAll(conSmirksReactants);
    				processed.addAll(conSmirksProducts);
    				
    			}
    			
    		}
			//apply bond and atom changes
    		applyChanges(aggregateProducts, changes,  aam, indexProduct);
    		//reset flags and get product separated
    		FragmentUtils.resetVisitedFlags(aggregateProducts);
    		IAtomContainerSet products = FragmentUtils.makeAtomContainerSet(aggregateProducts);
    		reaction.setReactants(reagents);
        	reaction.setProducts(products);

    	}
    	else {
    		System.err.println("Could not match the reaction");
    	}
		return reaction;
    	
    	//TODO check SMIRKS is correct (contains aam for all atoms)
    }

    
    /**
     * Apply a transform and get an unique product
     * @param reactantsSmiles
     * @param smirksString
     * @return
     * @throws CloneNotSupportedException
     * @throws CDKException
     */
    public IReaction applyTranform(String reactantsSmiles, String smirksString) throws CloneNotSupportedException, CDKException {
    	IReaction reaction = SilentChemObjectBuilder.getInstance().newInstance(IReaction.class);
    	SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
    	sp.kekulise(false);
    	IReaction smirksReaction = sp.parseReactionSmiles(smirksString);
    	String[] reactantsSmi = reactantsSmiles.split("\\.");
    	
    	IAtomContainerSet smirksReactants = smirksReaction.getReactants();
    	IAtomContainerSet smirksProducts = smirksReaction.getProducts();
    	IAtomContainer aggregateProducts = SilentChemObjectBuilder.getInstance().newAtomContainer();
    	IAtomContainerSet reagents = SilentChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
    	
    	//prepare products container to apply the reaction and  and index bonds
    	sp.kekulise(true);
    	int cpt = 1;
    	for (String reactant : reactantsSmi) {
    		IAtomContainer rea = sp.parseSmiles(reactant);
    		for (int i = 0; i < rea.getAtomCount(); i++) {
    			IAtom atom = rea.getAtom(i);
    			atom.setID(cpt+"");
    			atom.setProperty(CDKConstants.ATOM_ATOM_MAPPING, cpt);
    			cpt++;
    		}
    		aggregateProducts.add(rea);
    		reagents.addAtomContainer(rea.clone());
    	}

    	Map<IAtomContainer, int[]> mappings = new HashMap<IAtomContainer, int[]>();
    	//detect if all reactant patterns match
    	for (IAtomContainer smirksReactant : smirksReactants.atomContainers()) {
    		Pattern ptrn = Pattern.findSubstructure(smirksReactant);
    		int[] mapping = ptrn.match(aggregateProducts);
    		if (mapping.length > 0) {
    			mappings.put(smirksReactant, mapping);
    		}
    	}
    	
    	
    	//apply the reaction
    	if (!mappings.isEmpty()) {
    		//index atom product
    		Map<Integer,IAtom> indexProduct = new HashMap<Integer,IAtom>();
    		IAtomContainer aggregateSmirksProducts = SilentChemObjectBuilder.getInstance().newAtomContainer();
    		for (IAtomContainer smirksProduct : smirksProducts.atomContainers()) {
    			aggregateSmirksProducts.add(smirksProduct);
    			for (IAtom atom : smirksProduct.atoms()) {
    				indexProduct.put(atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING), atom);
    			}
    		}
    		
    		//apply changes
    		List<IBond> processed = new ArrayList<IBond>();
    		Map<IBond,IBond> changes = new HashMap<IBond,IBond>();
    		Map<Integer,IAtom> aam = new HashMap<Integer,IAtom>();
    		for (Entry<IAtomContainer, int[]> e : mappings.entrySet()) {
    			IAtomContainer smirksReactant = e.getKey();
    			int[] mapping = e.getValue();
    			for (int i = 0; i < mapping.length; i++) {
    				IAtom smirksReactantAtom = smirksReactant.getAtom(i);
    				IAtom smirksProductAtom = indexProduct.get(smirksReactantAtom.getProperty(CDKConstants.ATOM_ATOM_MAPPING));
    				aam.put(smirksReactantAtom.getProperty(CDKConstants.ATOM_ATOM_MAPPING), aggregateProducts.getAtom(mapping[i]));
    				
    				List<IBond> conSmirksReactants = smirksReactant.getConnectedBondsList(smirksReactantAtom);
    				List<IBond> conSmirksProducts = aggregateSmirksProducts.getConnectedBondsList(smirksProductAtom);
    				//remove process bond in reactants and products
    				for (IBond bond : new ArrayList<IBond>(conSmirksReactants)) {
    					if (processed.contains(bond)) {
    						conSmirksReactants.remove(bond);
    					}
    				}
    				for (IBond bond : new ArrayList<IBond>(conSmirksProducts)) {
    					if (processed.contains(bond)) {
    						conSmirksProducts.remove(bond);
    					}
    				}
    				
    				changes.putAll(findChanges(conSmirksReactants, conSmirksProducts));
    				
    				
    				processed.addAll(conSmirksReactants);
    				processed.addAll(conSmirksProducts);
    				
    			}
    			
    		}
			//apply bond and atom changes
    		applyChanges(aggregateProducts, changes,  aam, indexProduct);
    		//reset flags and get product separated
    		FragmentUtils.resetVisitedFlags(aggregateProducts);
    		IAtomContainerSet products = FragmentUtils.makeAtomContainerSet(aggregateProducts);
    		reaction.setReactants(reagents);
        	reaction.setProducts(products);

    	}
    	else {
    		System.err.println("Could not match the reaction");
    	}
    	
		return reaction;
    	
    	//TODO check SMIRKS is correct (contains aam for all atoms)
    }
    
    /**
     * Apply a transform and get all possible products
     * @param reagents
     * @param tranform
     * @return
     * @throws CloneNotSupportedException
     * @throws CDKException
     */
    public List<IReaction> applyTranform2(IAtomContainerSet reagents, IReaction tranform) throws CloneNotSupportedException, CDKException {    	
    	List<IReaction> results = new ArrayList<IReaction>();
    	List<IBitFingerprint> productFP = new ArrayList<IBitFingerprint>();
    	IAtomContainerSet smirksReactants = tranform.getReactants();
    	IAtomContainerSet smirksProducts = tranform.getProducts();
    	IAtomContainer aggregateProducts = SilentChemObjectBuilder.getInstance().newAtomContainer();
    	IAtomContainerSet reagentsCopy = SilentChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
    	
    	CircularFingerprinter fp = new CircularFingerprinter();
    	
    	//prepare products container to apply the reaction and  and index bonds
    	int cpt = 1;
    	for (IAtomContainer reactant : reagents.atomContainers()) {
    		reagentsCopy.addAtomContainer(reactant.clone());
    		aggregateProducts.add(reactant);
    		for (int i = 0; i < reactant.getAtomCount(); i++) {
    			IAtom atom = reactant.getAtom(i);
    			atom.setID(cpt+"");
    			atom.setProperty(CDKConstants.ATOM_ATOM_MAPPING, cpt);
    			cpt++;
    		}
    		cpt = 1;
    		for (int i = 0; i < reactant.getBondCount(); i++) {
    			IBond bond = reactant.getBond(i);
    			bond.setID(cpt+"");
    			bond.setProperty(CDKConstants.ATOM_ATOM_MAPPING, cpt);
    			cpt++;
    		}
    	}

    	Map<IAtomContainer, List<Map<IAtom, IAtom>> > mappings = new HashMap<IAtomContainer, List<Map<IAtom, IAtom>>>();
    	//detect if all reactant patterns match
    	for (IAtomContainer smirksReactant : smirksReactants.atomContainers()) {
    		Pattern ptrn = Pattern.findSubstructure(smirksReactant);
    		//query:target
    		List<Map<IAtom, IAtom>> mapping = getUniqueMappings(ptrn.matchAll(aggregateProducts).uniqueAtoms());
    		if (mapping.size() > 0) {
    			mappings.put(smirksReactant, mapping);
    		}
    	}
    	
    	
    	//apply the reaction
    	if (!mappings.isEmpty()) {
    		//index atom product
    		Map<Integer,IAtom> indexProduct = new HashMap<Integer,IAtom>();
    		IAtomContainer aggregateSmirksReactants = SilentChemObjectBuilder.getInstance().newAtomContainer();
    		IAtomContainer aggregateSmirksProducts = SilentChemObjectBuilder.getInstance().newAtomContainer();
    		for (IAtomContainer smirksReactant : smirksReactants.atomContainers()) {
    			aggregateSmirksReactants.add(smirksReactant);
    		}
    		for (IAtomContainer smirksProduct : smirksProducts.atomContainers()) {
    			aggregateSmirksProducts.add(smirksProduct);
    			for (IAtom atom : smirksProduct.atoms()) {
    				indexProduct.put(atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING), atom);
    			}
    		}
    		
    		
//    		store only the pattern which are matching multiple times
    		//Map<IAtomContainer,List<Map<IAtom, IAtom>>> multipleMatchPattern = new HashMap<IAtomContainer,List<Map<IAtom, IAtom>>>();

//    		store the index of the AAM in the List of solution corresponding to 1 pattern in order to generate all possible conditions
//    		ex: [P1={[M1,M2], P2=[M1,M2,M3]} ->allMultipleSolutionsMappings = [[M1,M2,M1,M2,M3]] -> List = [[0,1],[0,1,2]] 
//    		(list of lists to combine a mapping of one reactant to another one)
    		List<List<Integer>> indexOfTheMaps = new ArrayList<List<Integer>>();

//    		store the valid AAM (those which are unique)
    		List<Map<IAtom, IAtom>> validMappings = new ArrayList<Map<IAtom, IAtom>>();

//    		store all list of mappings which contains more than one solution in order to combine then and 
    		List<Map<IAtom, IAtom>> allMultipleSolutionsMappings = new ArrayList<Map<IAtom, IAtom>>();

    		int index = 0; 
    		for (Entry<IAtomContainer, List<Map<IAtom, IAtom>>> e1 : mappings.entrySet()) {
    			IAtomContainer pattern = e1.getKey();
    			List<Map<IAtom, IAtom>> mappingSolutions = e1.getValue();

    			if (mappingSolutions.size() == 1) {
    				validMappings.addAll(mappingSolutions);
    			}
    			else {
    				//multipleMatchPattern.put(pattern, mappingSolutions);
    				allMultipleSolutionsMappings.addAll(mappingSolutions);

    				
    				//add index of the map in a list and store it in indexOfTheMaps
    				List<Integer> list = new ArrayList<Integer>();
    				for (int i = 0; i < mappingSolutions.size(); i++) {
    					list.add(index);
    					index++;
    				}
    				indexOfTheMaps.add(list);
    			}
    		}
    		

    		//contains all mapping solutions which have to be applied to make the product
    		List<Map<IAtom, IAtom>> mappingSolutionToTest = new ArrayList<Map<IAtom, IAtom>>();

    		List<List<Integer>> mappingSolutionsToTest = new ArrayList<List<Integer>>();  
    		GenerateAllCombinationsFromMultipleLists(indexOfTheMaps, mappingSolutionsToTest, 0, new ArrayList<Integer>(), true);

    		//sort by size
    		Collections.sort(mappingSolutionsToTest, new Comparator<List<Integer>>(){
    			public int compare(List<Integer> a1, List<Integer> a2) {
    				return a1.size() - a2.size(); // smallest to biggest
    			}
    		});


    		for (List<Integer> solutionToTest : mappingSolutionsToTest) {
    			mappingSolutionToTest.clear();
    			mappingSolutionToTest.addAll(validMappings);

    			for (int elt : solutionToTest) {
    				mappingSolutionToTest.add(allMultipleSolutionsMappings.get(elt));
    			}

    	    	IReaction reaction = getInstance().newInstance(IReaction.class);
    			
    	    	//make reaction
        		List<IBond> processed = new ArrayList<IBond>();
        		Map<IBond,IBond> changes = new HashMap<IBond,IBond>();
        		Map<Integer,IAtom> aam = new HashMap<Integer,IAtom>();
        		for (Map<IAtom, IAtom> atomMap : mappingSolutionToTest) {
        			for (Entry<IAtom, IAtom> e : atomMap.entrySet()) {
        				IAtom smirksReactantAtom = e.getKey();
        				IAtom aggregateProductsAtom = e.getValue();
        				IAtom smirksProductAtom = indexProduct.get(smirksReactantAtom.getProperty(CDKConstants.ATOM_ATOM_MAPPING));

        				aam.put(smirksReactantAtom.getProperty(CDKConstants.ATOM_ATOM_MAPPING), aggregateProductsAtom);

        				List<IBond> conSmirksReactants = aggregateSmirksReactants.getConnectedBondsList(smirksReactantAtom);
        				List<IBond> conSmirksProducts = aggregateSmirksProducts.getConnectedBondsList(smirksProductAtom);
        				//remove process bond in reactants and products
        				for (IBond bond : new ArrayList<IBond>(conSmirksReactants)) {
        					if (processed.contains(bond)) {
        						conSmirksReactants.remove(bond);
        					}
        				}
        				for (IBond bond : new ArrayList<IBond>(conSmirksProducts)) {
        					if (processed.contains(bond)) {
        						conSmirksProducts.remove(bond);
        					}
        				}
        				
        				changes.putAll(findChanges(conSmirksReactants, conSmirksProducts));
        				
        				
        				processed.addAll(conSmirksReactants);
        				processed.addAll(conSmirksProducts);
        			}
        		}
    			IAtomContainer aggregateProductsCopy = aggregateProducts.clone();
        		//update aam
        		Map<Integer,IAtom> updatedAam = updateAam(aam, aggregateProductsCopy);
        		//uncomment to print bond changes
        		/*
        		for (Entry<IBond,IBond> ee : changes.entrySet()) {
        			IBond a = ee.getKey();
        			IBond b = ee.getValue();
        			if (a.getBegin() !=null && b.getBegin() != null)
        			System.out.println(a.getBegin().getSymbol()+":" +a.getOrder()+":" +a.getEnd().getSymbol() +"->" +
        					b.getBegin().getSymbol()+":" +b.getOrder()+":" +b.getEnd().getSymbol());
        			if (a.getBegin() !=null && b.getBegin() == null)
            			System.out.println(a.getBegin().getSymbol()+":" +a.getOrder()+":" +a.getEnd().getSymbol() +"->" +
            					"null");
        			if (a.getBegin() ==null && b.getBegin() != null)
        				System.out.println("null" +"->" +
            					b.getBegin().getSymbol()+":" +b.getOrder()+":" +b.getEnd().getSymbol());
        		}
        		*/
        		//apply bond and atom changes
        		try {
        			applyChanges(aggregateProductsCopy, changes,  updatedAam, indexProduct);
        		}
        		catch (Exception e) {
        			continue;
        		}
        		
        		// aggregateProductsCopy is null if a pattern match where a modification does not have to occur
        		// (ex: create a bond where a bond is already present)
        		if (!aggregateProductsCopy.isEmpty()) {
        			//Calculate FP and check if the product was already generated
            		fp.calculate(aggregateProductsCopy);
            		IBitFingerprint bfp = fp.getBitFingerprint(aggregateProductsCopy);
            		if (!productFP.contains(bfp)) {
    	        		//reset flags and get product separated
    	        		FragmentUtils.resetVisitedFlags(aggregateProductsCopy);
    	        		IAtomContainerSet products = FragmentUtils.makeAtomContainerSet(aggregateProductsCopy);
    	        		boolean valid = true;
    	        		//check if generated product is valid (only check modified atoms)
    	        		if (checkValence) {
    	        			for (IAtomContainer ac : products.atomContainers()) {
        	        			for (IAtom atom : updatedAam.values()) {
        	        				if (ac.contains(atom)) {
        	        					boolean isValencyValid = tools.isValencyOfMoleculeValid(ac, atom);
        	    	        			if (!isValencyValid){
        	        	        			valid = false;
        	        	        		}
        	        				}
        	        			}
        	        			//used previously to check all the atoms of all molecules
        	        			/*
        	        			boolean isValencyValid = tools.isValencyOfMoleculeValid(ac);
        	        			if (!isValencyValid){
            	        			valid = false;
            	        		}
            	        		*/
        	        		}
    	        		}
    	        		if (valid) {
    	        			reaction.setReactants(reagents);
        	            	reaction.setProducts(products);
        	            	results.add(reaction);
        	            	productFP.add(bfp);
    	        		}
            		}
        		}
    		}
    	}
    	else {
    		System.err.println("Could not match the reaction");
    	}
		return results;
    	
    	//TODO check SMIRKS is correct (contains aam for all atoms)
    }

    
    /**
     *  Apply a transform and get all possible products
     * @param reactantsSmiles
     * @param tranform
     * @return
     * @throws CloneNotSupportedException
     * @throws CDKException
     */
    public List<IReaction> applyTranform2(String reactantsSmiles, IReaction tranform) throws CloneNotSupportedException, CDKException {
    	List<IReaction> results = new ArrayList<IReaction>();
    	List<IBitFingerprint> productFP = new ArrayList<IBitFingerprint>();
    	SmilesParser sp = new SmilesParser(getInstance());
    	String[] reactantsSmi = reactantsSmiles.split("\\.");
    	
    	CircularFingerprinter fp = new CircularFingerprinter();

    	
    	IAtomContainerSet smirksReactants = tranform.getReactants();
    	IAtomContainerSet smirksProducts = tranform.getProducts();
    	IAtomContainer aggregateProducts = SilentChemObjectBuilder.getInstance().newAtomContainer();
    	IAtomContainerSet reagents = getInstance().newInstance(IAtomContainerSet.class);
    	
    	//prepare products container to apply the reaction and  and index bonds
    	sp.kekulise(true);
    	int cpt = 1;
    	int cpt2 = 1;
    	for (String reactant : reactantsSmi) {
    		IAtomContainer rea = sp.parseSmiles(reactant);
    		for (int i = 0; i < rea.getAtomCount(); i++) {
    			IAtom atom = rea.getAtom(i);
    			atom.setID(cpt+"");
    			atom.setProperty(CDKConstants.ATOM_ATOM_MAPPING, cpt);
    			cpt++;
    		}
    		for (int i = 0; i < rea.getBondCount(); i++) {
    			IBond bond = rea.getBond(i);
    			bond.setID(cpt2+"");
    			bond.setProperty(CDKConstants.ATOM_ATOM_MAPPING, cpt2);
    			cpt2++;
    		}
    		aggregateProducts.add(rea);
    		reagents.addAtomContainer(rea.clone());
    	}
    	
    	Map<IAtomContainer, List<Map<IAtom, IAtom>> > mappings = new HashMap<IAtomContainer, List<Map<IAtom, IAtom>>>();
    	//detect if all reactant patterns match
    	for (IAtomContainer smirksReactant : smirksReactants.atomContainers()) {
    		Pattern ptrn = Pattern.findSubstructure(smirksReactant);
    		//query:target
    		List<Map<IAtom, IAtom>> mapping = getUniqueMappings(ptrn.matchAll(aggregateProducts).uniqueAtoms());
    		if (mapping.size() > 0) {
    			mappings.put(smirksReactant, mapping);
    		}
    	}
    	
    	//apply the reaction
    	if (!mappings.isEmpty()) {
    		//index atom product
    		Map<Integer,IAtom> indexProduct = new HashMap<Integer,IAtom>();
    		IAtomContainer aggregateSmirksReactants = SilentChemObjectBuilder.getInstance().newAtomContainer();
    		IAtomContainer aggregateSmirksProducts = SilentChemObjectBuilder.getInstance().newAtomContainer();
    		for (IAtomContainer smirksReactant : smirksReactants.atomContainers()) {
    			aggregateSmirksReactants.add(smirksReactant);
    		}
    		for (IAtomContainer smirksProduct : smirksProducts.atomContainers()) {
    			aggregateSmirksProducts.add(smirksProduct);
    			for (IAtom atom : smirksProduct.atoms()) {
    				atom.setProperty("bcount", smirksProduct.getConnectedBondsCount(atom));
    				indexProduct.put(atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING), atom);
    			}
    		}
    		
    		
//    		store only the pattern which are matching multiple times
    		//Map<IAtomContainer,List<Map<IAtom, IAtom>>> multipleMatchPattern = new HashMap<IAtomContainer,List<Map<IAtom, IAtom>>>();

//    		store the index of the AAM in the List of solution corresponding to 1 pattern in order to generate all possible conditions
//    		ex: [P1={[M1,M2], P2=[M1,M2,M3]} ->allMultipleSolutionsMappings = [[M1,M2,M1,M2,M3]] -> List = [[0,1],[0,1,2]] 
//    		(list of lists to combine a mapping of one reactant to another one)
    		List<List<Integer>> indexOfTheMaps = new ArrayList<List<Integer>>();

//    		store the valid AAM (those which are unique)
    		List<Map<IAtom, IAtom>> validMappings = new ArrayList<Map<IAtom, IAtom>>();

//    		store all list of mappings which contains more than one solution in order to combine then and 
    		List<Map<IAtom, IAtom>> allMultipleSolutionsMappings = new ArrayList<Map<IAtom, IAtom>>();

    		int index = 0; 
    		for (Entry<IAtomContainer, List<Map<IAtom, IAtom>>> e1 : mappings.entrySet()) {
    			IAtomContainer pattern = e1.getKey();
    			List<Map<IAtom, IAtom>> mappingSolutions = e1.getValue();

    			if (mappingSolutions.size() == 1) {
//    				System.out.println("valid " + tools.makeSmiles(pattern, true, false, false));
    				validMappings.addAll(mappingSolutions);
    			}
    			else {
//    				System.out.println("unvalid " + tools.makeSmiles(pattern, true, false, false));
    				//multipleMatchPattern.put(pattern, mappingSolutions);
    				allMultipleSolutionsMappings.addAll(mappingSolutions);

    				
    				//add index of the map in a list and store it in indexOfTheMaps
    				List<Integer> list = new ArrayList<Integer>();
    				for (int i = 0; i < mappingSolutions.size(); i++) {
    					list.add(index);
    					index++;
    				}
    				indexOfTheMaps.add(list);
    			}
    		}
    		

    		//contains all mapping solutions which have to be applied to make the product
    		List<Map<IAtom, IAtom>> mappingSolutionToTest = new ArrayList<Map<IAtom, IAtom>>();

    		List<List<Integer>> mappingSolutionsToTest = new ArrayList<List<Integer>>();  
    		GenerateAllCombinationsFromMultipleLists(indexOfTheMaps, mappingSolutionsToTest, 0, new ArrayList<Integer>(), true);

    		//sort by size
    		Collections.sort(mappingSolutionsToTest, new Comparator<List<Integer>>(){
    			public int compare(List<Integer> a1, List<Integer> a2) {
    				return a1.size() - a2.size(); // smallest to biggest
    			}
    		});

    		for (List<Integer> solutionToTest : mappingSolutionsToTest) {
    			mappingSolutionToTest.clear();
    			mappingSolutionToTest.addAll(validMappings);

    			for (int elt : solutionToTest) {
    				mappingSolutionToTest.add(allMultipleSolutionsMappings.get(elt));
    			}
//    			System.out.println(solutionToTest + " " +mappingSolutionToTest.size() );
    	    	IReaction reaction = getInstance().newInstance(IReaction.class);
    			
    	    	//make reaction
        		
        		Map<IBond,IBond> changes = new HashMap<IBond,IBond>();
        		Map<Integer,List<IAtom>> aam = new HashMap<Integer,List<IAtom>>();
        		for (Map<IAtom, IAtom> atomMap : mappingSolutionToTest) {
        			//TODO MODIFIER AAM POUR AVOIR LIST 
//        			Set<IBond> processed = new HashSet<IBond>();
//        			System.out.println("-----------" + atomMap.size());
        			for (Entry<IAtom, IAtom> e : atomMap.entrySet()) {
//        				System.out.println( "size " + processed.size());
        				IAtom smirksReactantAtom = e.getKey();
        				IAtom aggregateProductsAtom = e.getValue();
        				
//        				System.out.println("id " + aggregateProductsAtom.getID());
        				IAtom smirksProductAtom = indexProduct.get(smirksReactantAtom.getProperty(CDKConstants.ATOM_ATOM_MAPPING));
        				
        				int mapVal = smirksReactantAtom.getProperty(CDKConstants.ATOM_ATOM_MAPPING);
        				if (aam.containsKey(mapVal)) {
        					List<IAtom> atoms = aam.get(mapVal);
        					atoms.add(aggregateProductsAtom);
        					aam.put(smirksReactantAtom.getProperty(CDKConstants.ATOM_ATOM_MAPPING), atoms);
        				}
        				else {
        					List<IAtom> atoms = new ArrayList<IAtom>();
        					atoms.add(aggregateProductsAtom);
        					aam.put(smirksReactantAtom.getProperty(CDKConstants.ATOM_ATOM_MAPPING), atoms);
        				}

        				List<IBond> conSmirksReactants = aggregateSmirksReactants.getConnectedBondsList(smirksReactantAtom);
        				List<IBond> conSmirksProducts = aggregateSmirksProducts.getConnectedBondsList(smirksProductAtom);
        				
        				//remove process bond in reactants and products
//        				for (IBond bond : new ArrayList<IBond>(conSmirksReactants)) {
//        					if (processed.contains(bond)) {
//        						conSmirksReactants.remove(bond);
//        					}
//        				}
//        				for (IBond bond : new ArrayList<IBond>(conSmirksProducts)) {
//        					if (processed.contains(bond)) {
//        						conSmirksProducts.remove(bond);
//        					}
//        				}
//        				
//        				processed.addAll(conSmirksReactants);
//        				processed.addAll(conSmirksProducts);
        				
        				Map<IBond,IBond> tempChanges = findChanges(conSmirksReactants, conSmirksProducts);
        				for (Entry<IBond,IBond> e2 : new HashMap<IBond,IBond>(tempChanges).entrySet()) {
        					IBond rBondSmirks = e2.getKey();
            	    		IBond pBondSmirks = e2.getValue();
            	    		if (pBondSmirks.getBegin() != null) {
            	    			IAtom begin = pBondSmirks.getBegin();
            	    			IAtom end = pBondSmirks.getEnd();
            	    			if (!aam.containsKey(begin.getProperty(CDKConstants.ATOM_ATOM_MAPPING)) ||
                	    				!aam.containsKey(end.getProperty(CDKConstants.ATOM_ATOM_MAPPING))) {
//            	    				System.out.println("ps " + 
//                	    				begin.getProperty(CDKConstants.ATOM_ATOM_MAPPING) + " " +
//                	    				end.getProperty(CDKConstants.ATOM_ATOM_MAPPING)
//                	    						+aam.keySet());
            	    				tempChanges.remove(rBondSmirks);
                	    		}
            	    			else {
            	    				pBondSmirks.setProperty("ref", smirksProductAtom);
//            	    				if (smirksProductAtom.equals(begin)) {
//            	    					
//            	    				}
//            	    				else {
//            	    					
//            	    				}
            	    			}
            	    		}
        				}
        				changes.putAll(tempChanges);
        				
//        				System.out.println(smirksReactantAtom.getSymbol() + " " +conSmirksReactants.size() + " " + processed.size() + " " + changes.size());
//        				System.out.println(smirksReactantAtom.getSymbol() + " " +conSmirksReactants.size() + " " + changes.size());

//        				System.out.println(conSmirksProducts);
        				
        			}
        		}
    			IAtomContainer aggregateProductsCopy = aggregateProducts.clone();
        		//update aam
        		Map<Integer,List<IAtom>> updatedAam = updateAam2(aam, aggregateProductsCopy);
        		//uncomment to print bond changes
        		/*
        		//System.out.println("CHANGED " + changes.size());
        		for (Entry<IBond,IBond> ee : changes.entrySet()) {
        			IBond a = ee.getKey();
        			IBond b = ee.getValue();
        			if (a.getBegin() !=null && b.getBegin() != null)
        			System.out.println(a.getBegin().getSymbol()+":" +a.getOrder()+":" +a.getEnd().getSymbol() +"->" +
        					b.getBegin().getSymbol()+":" +b.getOrder()+":" +b.getEnd().getSymbol());
        			if (a.getBegin() !=null && b.getBegin() == null)
            			System.out.println(a.getBegin().getSymbol()+":" +a.getOrder()+":" +a.getEnd().getSymbol() +"->" +
            					"null");
        			if (a.getBegin() ==null && b.getBegin() != null)
        				System.out.println("null" +"->" +
            					b.getBegin().getSymbol()+":" +b.getOrder()+":" +b.getEnd().getSymbol());
        		}
        		*/
        		
        		//apply bond and atom changes
        		//try {
        			applyChanges2(aggregateProductsCopy, changes,  updatedAam, indexProduct);
//        		}
//        		catch (Exception e) {
//        			continue;
//        		}
//        		// aggregateProductsCopy is null if a pattern match where a modification does not have to occur
        		// (ex: create a bond where a bond is already present)
        		if (!aggregateProductsCopy.isEmpty()) {
        			//Calculate FP and check if the product was already generated
            		fp.calculate(aggregateProductsCopy);
            		IBitFingerprint bfp = fp.getBitFingerprint(aggregateProductsCopy);
            		if (!productFP.contains(bfp)) {
    	        		//reset flags and get product separated
    	        		FragmentUtils.resetVisitedFlags(aggregateProductsCopy);
    	        		IAtomContainerSet products = FragmentUtils.makeAtomContainerSet(aggregateProductsCopy);
    	        		boolean valid = true;
    	        		//check if generated product is valid (only check modified atoms)
    	        		if (checkValence) {
    	        			for (IAtomContainer ac : products.atomContainers()) {
        	        			for (List<IAtom> atoms : updatedAam.values()) {
        	        				for (IAtom atom : atoms) {
        	        					if (ac.contains(atom)) {
            	        					boolean isValencyValid = tools.isValencyOfMoleculeValid(ac, atom);
            	    	        			if (!isValencyValid){
            	        	        			valid = false;
            	        	        		}
            	        				}
        	        				}
        	        				
        	        			}
        	        			//used previously to check all the atoms of all molecules
        	        			/*
        	        			boolean isValencyValid = tools.isValencyOfMoleculeValid(ac);
        	        			if (!isValencyValid){
            	        			valid = false;
            	        		}
            	        		*/
        	        		}
    	        		}
    	        		if (valid) {
    	        			reaction.setReactants(reagents);
        	            	reaction.setProducts(products);
        	            	results.add(reaction);
        	            	productFP.add(bfp);
    	        		}
            		}
        		}
    		}
    	}
    	else {
    		System.err.println("Could not match the reaction");
    	}
		return results;
    	
    	//TODO check SMIRKS is correct (contains aam for all atoms)
    }

    /**
     *  Apply a transform and get all possible products
     * @param reactantsSmiles
     * @param smirksString
     * @return
     * @throws CloneNotSupportedException
     * @throws CDKException
     */
    public List<IReaction> applyTranform2(String reactantsSmiles, String smirksString) throws CloneNotSupportedException, CDKException {
    	List<IReaction> results = new ArrayList<IReaction>();
    	List<IBitFingerprint> productFP = new ArrayList<IBitFingerprint>();
    	SmilesParser sp = new SmilesParser(getInstance());
    	sp.kekulise(false);
    	IReaction smirksReaction = sp.parseReactionSmiles(smirksString);
    	String[] reactantsSmi = reactantsSmiles.split("\\.");
    	
    	CircularFingerprinter fp = new CircularFingerprinter();
    	
    	IAtomContainerSet smirksReactants = smirksReaction.getReactants();
    	IAtomContainerSet smirksProducts = smirksReaction.getProducts();
    	IAtomContainer aggregateProducts = SilentChemObjectBuilder.getInstance().newAtomContainer();
    	IAtomContainerSet reagents = getInstance().newInstance(IAtomContainerSet.class);
    	
    	//prepare products container to apply the reaction and  and index bonds
    	sp.kekulise(true);
    	int cpt = 1;
    	for (String reactant : reactantsSmi) {
    		IAtomContainer rea = sp.parseSmiles(reactant);
    		for (int i = 0; i < rea.getAtomCount(); i++) {
    			IAtom atom = rea.getAtom(i);
    			atom.setID(cpt+"");
    			atom.setProperty(CDKConstants.ATOM_ATOM_MAPPING, cpt);
    			cpt++;
    		}
    		cpt = 1;
    		for (int i = 0; i < rea.getBondCount(); i++) {
    			IBond bond = rea.getBond(i);
    			bond.setID(cpt+"");
    			bond.setProperty(CDKConstants.ATOM_ATOM_MAPPING, cpt);
    			cpt++;
    		}
    		aggregateProducts.add(rea);
    		reagents.addAtomContainer(rea.clone());
    	}

    	Map<IAtomContainer, List<Map<IAtom, IAtom>> > mappings = new HashMap<IAtomContainer, List<Map<IAtom, IAtom>>>();
    	//detect if all reactant patterns match
    	for (IAtomContainer smirksReactant : smirksReactants.atomContainers()) {
    		Pattern ptrn = Pattern.findSubstructure(smirksReactant);
    		//query:target
    		List<Map<IAtom, IAtom>> mapping = getUniqueMappings(ptrn.matchAll(aggregateProducts).uniqueAtoms());
    		if (mapping.size() > 0) {
    			mappings.put(smirksReactant, mapping);
    		}
    	}
    	
    	
    	//apply the reaction
    	if (!mappings.isEmpty()) {
    		//index atom product
    		Map<Integer,IAtom> indexProduct = new HashMap<Integer,IAtom>();
    		IAtomContainer aggregateSmirksReactants = SilentChemObjectBuilder.getInstance().newAtomContainer();
    		IAtomContainer aggregateSmirksProducts = SilentChemObjectBuilder.getInstance().newAtomContainer();
    		for (IAtomContainer smirksReactant : smirksReactants.atomContainers()) {
    			aggregateSmirksReactants.add(smirksReactant);
    		}
    		for (IAtomContainer smirksProduct : smirksProducts.atomContainers()) {
    			aggregateSmirksProducts.add(smirksProduct);
    			for (IAtom atom : smirksProduct.atoms()) {
    				indexProduct.put(atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING), atom);
    			}
    		}
    		
    		
//    		store only the pattern which are matching multiple times
    		//Map<IAtomContainer,List<Map<IAtom, IAtom>>> multipleMatchPattern = new HashMap<IAtomContainer,List<Map<IAtom, IAtom>>>();

//    		store the index of the AAM in the List of solution corresponding to 1 pattern in order to generate all possible conditions
//    		ex: [P1={[M1,M2], P2=[M1,M2,M3]} ->allMultipleSolutionsMappings = [[M1,M2,M1,M2,M3]] -> List = [[0,1],[0,1,2]] 
//    		(list of lists to combine a mapping of one reactant to another one)
    		List<List<Integer>> indexOfTheMaps = new ArrayList<List<Integer>>();

//    		store the valid AAM (those which are unique)
    		List<Map<IAtom, IAtom>> validMappings = new ArrayList<Map<IAtom, IAtom>>();

//    		store all list of mappings which contains more than one solution in order to combine then and 
    		List<Map<IAtom, IAtom>> allMultipleSolutionsMappings = new ArrayList<Map<IAtom, IAtom>>();

    		int index = 0; 
    		for (Entry<IAtomContainer, List<Map<IAtom, IAtom>>> e1 : mappings.entrySet()) {
    			IAtomContainer pattern = e1.getKey();
    			List<Map<IAtom, IAtom>> mappingSolutions = e1.getValue();

    			if (mappingSolutions.size() == 1) {
    				validMappings.addAll(mappingSolutions);
    			}
    			else {
    				//multipleMatchPattern.put(pattern, mappingSolutions);
    				allMultipleSolutionsMappings.addAll(mappingSolutions);

    				
    				//add index of the map in a list and store it in indexOfTheMaps
    				List<Integer> list = new ArrayList<Integer>();
    				for (int i = 0; i < mappingSolutions.size(); i++) {
    					list.add(index);
    					index++;
    				}
    				indexOfTheMaps.add(list);
    			}
    		}
    		

    		//contains all mapping solutions which have to be applied to make the product
    		List<Map<IAtom, IAtom>> mappingSolutionToTest = new ArrayList<Map<IAtom, IAtom>>();

    		List<List<Integer>> mappingSolutionsToTest = new ArrayList<List<Integer>>();  
    		GenerateAllCombinationsFromMultipleLists(indexOfTheMaps, mappingSolutionsToTest, 0, new ArrayList<Integer>(), true);

    		//sort by size
    		Collections.sort(mappingSolutionsToTest, new Comparator<List<Integer>>(){
    			public int compare(List<Integer> a1, List<Integer> a2) {
    				return a1.size() - a2.size(); // smallest to biggest
    			}
    		});


    		for (List<Integer> solutionToTest : mappingSolutionsToTest) {
    			mappingSolutionToTest.clear();
    			mappingSolutionToTest.addAll(validMappings);

    			for (int elt : solutionToTest) {
    				mappingSolutionToTest.add(allMultipleSolutionsMappings.get(elt));
    			}

    	    	IReaction reaction = getInstance().newInstance(IReaction.class);
    			
    	    	//make reaction
        		List<IBond> processed = new ArrayList<IBond>();
        		Map<IBond,IBond> changes = new HashMap<IBond,IBond>();
        		Map<Integer,IAtom> aam = new HashMap<Integer,IAtom>();
        		for (Map<IAtom, IAtom> atomMap : mappingSolutionToTest) {
        			for (Entry<IAtom, IAtom> e : atomMap.entrySet()) {
        				IAtom smirksReactantAtom = e.getKey();
        				IAtom aggregateProductsAtom = e.getValue();
        				IAtom smirksProductAtom = indexProduct.get(smirksReactantAtom.getProperty(CDKConstants.ATOM_ATOM_MAPPING));

        				aam.put(smirksReactantAtom.getProperty(CDKConstants.ATOM_ATOM_MAPPING), aggregateProductsAtom);

        				List<IBond> conSmirksReactants = aggregateSmirksReactants.getConnectedBondsList(smirksReactantAtom);
        				List<IBond> conSmirksProducts = aggregateSmirksProducts.getConnectedBondsList(smirksProductAtom);
        				//remove process bond in reactants and products
        				for (IBond bond : new ArrayList<IBond>(conSmirksReactants)) {
        					if (processed.contains(bond)) {
        						conSmirksReactants.remove(bond);
        					}
        				}
        				for (IBond bond : new ArrayList<IBond>(conSmirksProducts)) {
        					if (processed.contains(bond)) {
        						conSmirksProducts.remove(bond);
        					}
        				}
        				
        				changes.putAll(findChanges(conSmirksReactants, conSmirksProducts));
        				
        				
        				processed.addAll(conSmirksReactants);
        				processed.addAll(conSmirksProducts);
        			}
        		}
    			IAtomContainer aggregateProductsCopy = aggregateProducts.clone();
        		//update aam
        		Map<Integer,IAtom> updatedAam = updateAam(aam, aggregateProductsCopy);
        		//uncomment to print bond changes
        		/*
        		for (Entry<IBond,IBond> ee : changes.entrySet()) {
        			IBond a = ee.getKey();
        			IBond b = ee.getValue();
        			if (a.getBegin() !=null && b.getBegin() != null)
        			System.out.println(a.getBegin().getSymbol()+":" +a.getOrder()+":" +a.getEnd().getSymbol() +"->" +
        					b.getBegin().getSymbol()+":" +b.getOrder()+":" +b.getEnd().getSymbol());
        			if (a.getBegin() !=null && b.getBegin() == null)
            			System.out.println(a.getBegin().getSymbol()+":" +a.getOrder()+":" +a.getEnd().getSymbol() +"->" +
            					"null");
        			if (a.getBegin() ==null && b.getBegin() != null)
        				System.out.println("null" +"->" +
            					b.getBegin().getSymbol()+":" +b.getOrder()+":" +b.getEnd().getSymbol());
        		}
        		*/
        		//apply bond and atom changes
        		try {
        			applyChanges(aggregateProductsCopy, changes,  updatedAam, indexProduct);
        		}
        		catch (Exception e) {
        			continue;
        		}

        		// aggregateProductsCopy is null if a pattern match where a modification does not have to occur
        		// (ex: create a bond where a bond is already present)
        		if (!aggregateProductsCopy.isEmpty()) {
        			//Calculate FP and check if the product was already generated
            		fp.calculate(aggregateProductsCopy);
            		IBitFingerprint bfp = fp.getBitFingerprint(aggregateProductsCopy);
            		if (!productFP.contains(bfp)) {
    	        		//reset flags and get product separated
    	        		FragmentUtils.resetVisitedFlags(aggregateProductsCopy);
    	        		IAtomContainerSet products = FragmentUtils.makeAtomContainerSet(aggregateProductsCopy);
    	        		boolean valid = true;
    	        		//check if generated product is valid (only check modified atoms)
    	        		if (checkValence) {
    	        			for (IAtomContainer ac : products.atomContainers()) {
        	        			for (IAtom atom : updatedAam.values()) {
        	        				if (ac.contains(atom)) {
        	        					boolean isValencyValid = tools.isValencyOfMoleculeValid(ac, atom);
        	    	        			if (!isValencyValid){
        	        	        			valid = false;
        	        	        		}
        	        				}
        	        			}
        	        			//used previously to check all the atoms of all molecules
        	        			/*
        	        			boolean isValencyValid = tools.isValencyOfMoleculeValid(ac);
        	        			if (!isValencyValid){
            	        			valid = false;
            	        		}
            	        		*/
        	        		}
    	        		}
    	        		if (valid) {
    	        			reaction.setReactants(reagents);
        	            	reaction.setProducts(products);
        	            	results.add(reaction);
        	            	productFP.add(bfp);
    	        		}
            		}
        		}
    		}
    	}
    	else {
    		System.err.println("Could not match the reaction");
    	}
		return results;
    	
    	//TODO check SMIRKS is correct (contains aam for all atoms)
    }

    
	/**
	 * Get each unique Mapping Solution, i.e. get the different AAM combination
	 * @param m
	 * @return
	 */
	private List<Map<IAtom,IAtom>> getUniqueMappings(Mappings m) {
		List<Map<IAtom,IAtom>> res = new ArrayList<Map<IAtom,IAtom>>();
//		List<List<IAtom>> exclusiveMapping = new ArrayList<List<IAtom>>();
		List<List<String>> uniqueMapping = new ArrayList<List<String>>();
		for (Map<IAtom,IAtom> map : m.toAtomMap()) {
			List<String> mappedAtom = new ArrayList<String>();
			for (Entry<IAtom,IAtom> e : map.entrySet()) {
				mappedAtom.add(e.getValue().getID());
			}
			Collections.sort(mappedAtom);
			if (!uniqueMapping.contains(mappedAtom)) {
				uniqueMapping.add(mappedAtom);
				res.add(map);
			}
				
		}
		return res;
	}
    
	/**
	 * Get each unique Mapping Solution, i.e. get the different AAM combination
	 * @param m
	 * @return
	 */
	private List<Map<IBond,IBond>> getUniqueBondMappings(Mappings m) {
		List<Map<IBond,IBond>> res = new ArrayList<Map<IBond,IBond>>();
//		List<List<IAtom>> exclusiveMapping = new ArrayList<List<IAtom>>();
		List<List<String>> uniqueMapping = new ArrayList<List<String>>();
		for (Map<IBond,IBond> map : m.toBondMap()) {
			List<String> mappedAtom = new ArrayList<String>();
			for (Entry<IBond,IBond> e : map.entrySet()) {
				mappedAtom.add(e.getValue().getID());
			}
			Collections.sort(mappedAtom);
			if (!uniqueMapping.contains(mappedAtom)) {
				uniqueMapping.add(mappedAtom);
				res.add(map);
			}
				
		}
		return res;
	}
	
	/**
	 * Generate combination from multiple lists
	 * X: [1, 2] 
	 * Y: [3, 4, 5, 6]
	 * Result: [[1, 3], [1, 4], [1, 5], [1, 6], [2, 3], [2, 4], [2, 5], [2, 6]]
	 * if parameter intraCombinations is true, generate combination between the lists and for each combination 
	 * generate all permutations inside the permuted list
	 * X: [1, 2] 
	 * [3, 4, 5, 6]
	 * Result: [[1, 3], [1, 4], [1, 5], [1, 6], [2, 3], [2, 4], [2, 5], [2, 6], [1, 4, 3], [1, 5, 3], [1, 5, 4], [1, 6, 3], 
	 * [1, 6, 4], [1, 6, 5], [2, 4, 3], [2, 5, 3], [2, 5, 4], [2, 6, 3], [2, 6, 4], [2, 6, 5], [2, 1, 3], [2, 1, 4], [2, 1, 5], 
	 * [2, 1, 6], [1, 5, 3, 4], [1, 6, 3, 4], [1, 6, 3, 5], [1, 6, 4, 5], [2, 5, 3, 4], [2, 6, 3, 4], [2, 6, 3, 5], [2, 6, 4, 5], 
	 * [2, 1, 4, 3], [2, 1, 5, 3], [2, 1, 5, 4], [2, 1, 6, 3], [2, 1, 6, 4], [2, 1, 6, 5], [1, 6, 3, 4, 5], [2, 6, 3, 4, 5], 
	 * [2, 1, 5, 3, 4], [2, 1, 6, 3, 4], [2, 1, 6, 3, 5], [2, 1, 6, 4, 5], [2, 1, 6, 3, 4, 5]]
	 * https://stackoverflow.com/questions/17192796/generate-all-combinations-from-multiple-lists
	 * https://javahungry.blogspot.com/2014/02/algorithm-for-combinations-of-string-java-code-with-example.html
	 * https://www.geeksforgeeks.org/print-all-possible-combinations-of-r-elements-in-a-given-array-of-size-n/
	 * @param Lists
	 * @param result
	 * @param depth
	 * @param current
	 */
	private void GenerateAllCombinationsFromMultipleLists(List<List<Integer>> Lists, List<List<Integer>> result, int depth, 
			List<Integer> current, boolean intraCombinations) {
		if(depth == Lists.size()) {
			result.add(new ArrayList<Integer>(current));
			//System.out.println("res " + result);
			return;
		}

		for(int i = 0; i < Lists.get(depth).size(); i++) {
			int t = Lists.get(depth).get(i);

			List<Integer> temp = new ArrayList<Integer>(current);
			temp.add(t);

			GenerateAllCombinationsFromMultipleLists(Lists, result, depth + 1, temp, intraCombinations);

			//generate all combination inside a list and return a list which contains lists of different combination
			//ex sublist [3, 4, 5] -> combinationResultsOfTheSublist [[3, 4, 5], [3, 4], [3, 5], [3], [4, 5], [4], [5]]
			//limit to 5 times (stoichiometry)
			if (intraCombinations == true && i > 0 && i < stoichiometry) {
				List<Integer> sublist = Lists.get(depth).subList(0, i); 
				Integer[] sublistArray = sublist.toArray(new Integer[sublist.size()]);
				List<List<Integer>> combinationResultsOfTheSublist = new ArrayList<List<Integer>>();
				generateAllCombinationsInsideAList(sublistArray, combinationResultsOfTheSublist);

				//generate all combinations of the previously generated results with the other list in Lists
				for(int j = 0; j < combinationResultsOfTheSublist.size(); j++) {
					List<Integer> temp2 = new ArrayList<Integer>(temp);
					temp2.addAll(combinationResultsOfTheSublist.get(j));
					GenerateAllCombinationsFromMultipleLists(Lists, result, depth + 1, temp2, intraCombinations);
				}
			}

		}
	}
	
	/**
	 * @param elems
	 * @param result
	 */
	private void generateAllCombinationsInsideAList(Integer[] elems, List<List<Integer>> result) {
	    List<Integer> combination = new ArrayList<>();
	    recursive_combinations(combination, 0, elems, result);
	}
	
	/**
	 * @param combination
	 * @param ndx
	 * @param elems
	 * @param result
	 */
	private void recursive_combinations(List<Integer> combination,
                                      int ndx, Integer[] elems, List<List<Integer>> result) {
	    if(ndx == elems.length) {
	
	        // (reached end of list after selecting/not selecting)
	        if (!combination.isEmpty())
	            result.add(new ArrayList<Integer>(combination));
	
	    } 
	    else {
	
	        // (include element at ndx)
	        combination.add(elems[ndx]);
	        recursive_combinations(combination, ndx+1, elems, result);
	
	        // (don't include element at ndx)
	        combination.remove(elems[ndx]);
	        recursive_combinations(combination, ndx+1, elems, result);
	
	    }
	}
	
	private String makeBondID(IAtom begin, IAtom end) {
		int aam1 = (int) begin.getProperty(CDKConstants.ATOM_ATOM_MAPPING);
		int aam2 = (int) end.getProperty(CDKConstants.ATOM_ATOM_MAPPING);
		if  (aam1 < aam2) {
			return aam1 + "-" + aam2;
		}
		else {
			return aam2 + "-" + aam1;
		}
	}
	private Map<Integer,IAtom> updateAam(Map<Integer,IAtom> aam, IAtomContainer ac) {
		Map<Integer,IAtom> newAam = new HashMap<Integer,IAtom>(aam);
		Map<Integer,IAtom> temp = new HashMap<Integer,IAtom>();
		for (int i = 0; i < ac.getAtomCount(); i++) {
			IAtom atom = ac.getAtom(i);
			temp.put(atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING), atom);
		}
		for (Entry<Integer,IAtom> e : aam.entrySet()) {
			int key = e.getKey();
			IAtom value = e.getValue();
			newAam.put(key, temp.get(value.getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
		}
		return newAam;
	}
	
	private Map<Integer,List<IAtom>> updateAam2(Map<Integer,List<IAtom>> aam, IAtomContainer ac) {
		Map<Integer,List<IAtom>> newAam = new HashMap<Integer,List<IAtom>>();
		Map<Integer,IAtom> temp = new HashMap<Integer,IAtom>();
		for (int i = 0; i < ac.getAtomCount(); i++) {
			IAtom atom = ac.getAtom(i);
			temp.put(atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING), atom);
		}
		for (Entry<Integer,List<IAtom>> e : aam.entrySet()) {
			int key = e.getKey();
			List<IAtom> values  = new ArrayList<IAtom>(e.getValue());
			for (IAtom atom : values) {
//				System.out.println(key + " " + atom.getSymbol() + " " + atom.getProperties());
				int mapVal = atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING);
				if (newAam.containsKey(key)) {
					List<IAtom> atoms = newAam.get(key);
					atoms.add(temp.get(mapVal));
					newAam.put(key, atoms);
				}
				else {
					List<IAtom> atoms = new ArrayList<IAtom>();
					atoms.add(temp.get(mapVal));
					newAam.put(key, atoms);
				}
			}
		}
		return newAam;
	}
	
	private Map<String,List<IBond>> updateBbm(Map<String,List<IBond>> bbm, IAtomContainer ac) {
		Map<String,List<IBond>> newBbm = new HashMap<String,List<IBond>>();
		Map<Integer,IBond> temp = new HashMap<Integer,IBond>();
		for (int i = 0; i < ac.getBondCount(); i++) {
			IBond bond = ac.getBond(i);
			temp.put(bond.getProperty(CDKConstants.ATOM_ATOM_MAPPING), bond);
		}
		for (Entry<String,List<IBond>> e : bbm.entrySet()) {
			String key = e.getKey();
			List<IBond> values  = new ArrayList<IBond>(e.getValue());
			for (IBond atom : values) {
				int mapVal = atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING);
				if (newBbm.containsKey(key)) {
					List<IBond> atoms = newBbm.get(key);
					atoms.add(temp.get(mapVal));
					newBbm.put(key, atoms);
				}
				else {
					List<IBond> atoms = new ArrayList<IBond>();
					atoms.add(temp.get(mapVal));
					newBbm.put(key, atoms);
				}
			}
		}
		return newBbm;
	}
	
    private Map<IBond,IBond> findChanges(List<IBond> inReactants, List<IBond> inProducts) {
    	Map<IBond,IBond> res = new HashMap<IBond,IBond>();
    	List<IBond> rMatched = new ArrayList<IBond>();
    	List<IBond> pMatched = new ArrayList<IBond>();
    	for (IBond rBond : inReactants) {
    		String rid = "";
    		if ((int)rBond.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING) > 
    			(int)rBond.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING)) {
    			rid = rBond.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING) + "-" + 
    					rBond.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING);
    		}
    		else {
    			rid = rBond.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING) + "-" + 
    					rBond.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING);
    		}
    		for (IBond pBond : inProducts) {
    			if (pMatched.contains(pBond))
    				continue;
    			String pid = "";
        		if ((int)pBond.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING) > 
        			(int)pBond.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING)) {
        			pid = pBond.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING) + "-" + 
        					pBond.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING);
        		}
        		else {
        			pid = pBond.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING) + "-" + 
        					pBond.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING);
        		}
        		if (rid.equals(pid)) {
        			rMatched.add(rBond);
        			pMatched.add(pBond);
        			if (!areBothBondsEquals(rBond, pBond)) {
        				res.put(rBond, pBond);
        			}
        		}
    		}
    	}
    	
    	//made bond
    	List<IBond> formed = inProducts;
    	formed.removeAll(pMatched);
    	for (IBond bond : formed) {
    		res.put(new Bond(), bond);
    	}
    	
    	//broken bond
    	List<IBond> broken = inReactants;
    	broken.removeAll(rMatched);
    	for (IBond bond : broken) {
    		res.put(bond, new Bond());
    	}
		return res;
    }
    
    private boolean areBothBondsEquals(IBond rBond, IBond pBond) {
    	boolean isEqual = true;
    	if (rBond.isAromatic() != pBond.isAromatic()) {
    		//rBond.setIsAromatic(pBond.isAromatic());
    		isEqual = false;
    	}
    	else {
    		if (!rBond.getOrder().equals(pBond.getOrder())) {
    			//rBond.setOrder(pBond.getOrder());
    			isEqual = false;
    		}
    	}
    	if (!rBond.getStereo().equals(pBond.getStereo())) {
    		//rBond.setStereo(pBond.getStereo());
    		isEqual = false;
    	}
    	
    	IAtom pBeginAtom;
    	IAtom pEndAtom;
    	if (rBond.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING).equals(pBond.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING))) {
    		pBeginAtom = pBond.getBegin();
    		pEndAtom = pBond.getEnd();
    	}
    	else {
    		pBeginAtom = pBond.getEnd();
    		pEndAtom = pBond.getBegin();
    	}
    	
    	if (rBond.getBegin().getFormalCharge() != pBeginAtom.getFormalCharge()) {
    		//rBond.getBegin().setFormalCharge(pBeginAtom.getFormalCharge());
    		isEqual = false;
    	}
    	if (rBond.getEnd().getFormalCharge() != pEndAtom.getFormalCharge()) {
    		//rBond.getEnd().setFormalCharge(pEndAtom.getFormalCharge());
    		isEqual = false;
    	}
    	if (rBond.getBegin().getCharge() != pBeginAtom.getCharge()) {
    		//rBond.getBegin().setCharge(pBeginAtom.getCharge());
    		isEqual = false;
    	}
    	if (rBond.getEnd().getCharge() != pEndAtom.getCharge()) {
    		//rBond.getEnd().setCharge(pEndAtom.getCharge());
    		isEqual = false;
    	}
    	if (rBond.getBegin().getExactMass() != pBeginAtom.getExactMass()) {
    		//rBond.getBegin().setExactMass(pBeginAtom.getExactMass());
    		isEqual = false;
    	}
    	if (rBond.getEnd().getExactMass() != pEndAtom.getExactMass()) {
    		//rBond.getEnd().setExactMass(pEndAtom.getExactMass());
    		isEqual = false;
    	}
    	if (rBond.getBegin().isInRing() != pBeginAtom.isInRing()) {
    		//rBond.getBegin().setIsInRing(pBeginAtom.isInRing());
    		isEqual = false;
    	}
    	if (rBond.getEnd().isInRing() != pEndAtom.isInRing()) {
    		//rBond.getEnd().setIsInRing(pEndAtom.isInRing());
    		isEqual = false;
    	}
//    	if (rBond.getBegin().getImplicitHydrogenCount() != pBeginAtom.getImplicitHydrogenCount()) {
//    		rBond.getBegin().setImplicitHydrogenCount(pBeginAtom.getImplicitHydrogenCount());
//    		isEqual = false;
//    	}
//    	if (rBond.getEnd().getImplicitHydrogenCount() != pEndAtom.getImplicitHydrogenCount()) {
//    		rBond.getEnd().setImplicitHydrogenCount(pEndAtom.getImplicitHydrogenCount());
//    		isEqual = false;
//    	}
//
//    	if (rBond.getBegin().getHybridization() != pBeginAtom.getHybridization()) {
//    		rBond.getBegin().setHybridization(pBeginAtom.getHybridization());
//    		isEqual = false;
//    	}
//    	if (rBond.getEnd().getHybridization() != pEndAtom.getHybridization()) {
//    		rBond.getEnd().setHybridization(pEndAtom.getHybridization());
//    		isEqual = false;
//    	}
//    	if (rBond.getBegin().getValency()!= pBeginAtom.getValency()) {
//    		rBond.getBegin().setValency(pBeginAtom.getValency());
//    		isEqual = false;
//    	}
//    	if (rBond.getEnd().getValency() != pEndAtom.getValency()) {
//    		rBond.getEnd().setValency(pEndAtom.getValency());
//    		isEqual = false;
//    	}
    	return isEqual;
    }
    
    private void applyChanges(IAtomContainer ac, Map<IBond,IBond> changesMap, 
    		Map<Integer,IAtom> aam, Map<Integer,IAtom> indexProduct) {
    	Set<String> processed = new HashSet<String>();
    	//si deja process ne pas ajouter ou modifier
    	for (Entry<IBond,IBond> e : changesMap.entrySet()) {
    		IBond rBondSmirks = e.getKey();
    		IBond pBondSmirks = e.getValue();
    		IAtom begin;
    		IAtom end;
    		//made
    		if (rBondSmirks.getBegin() == null) {
    			IBond newBond = new Bond();
    			begin = aam.get(pBondSmirks.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING));
    			end = aam.get(pBondSmirks.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING));
    			if (ac.getBond(begin, end) != null) {
    				ac.removeAllElements();
    				return;
    			}
    			newBond.setAtom(begin, 0);
    			newBond.setAtom(end, 1);
    			newBond.setIsAromatic(pBondSmirks.isAromatic());
    			newBond.setOrder(pBondSmirks.getOrder());
    			newBond.setStereo(pBondSmirks.getStereo());
    			ac.addBond(newBond);
    			applyAtomModification(begin, indexProduct.get(pBondSmirks.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
    			applyAtomModification(end, indexProduct.get(pBondSmirks.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
    			
    			//correct hydrogen
    			tools.caclulateHydrogen(ac, begin, true);
    			tools.caclulateHydrogen(ac, end, true);
    			//OLD Method
    			/*
    			int hBegin = tools.caclulateHydrogen2(ac, begin, true);
    			if (hBegin > -1) 
    				begin.setImplicitHydrogenCount(hBegin);
    			int hEnd = tools.caclulateHydrogen2(ac, end, true);
    			if (hEnd > -1) 
    				end.setImplicitHydrogenCount(hEnd)
    			*/
    		}
    		//broken
    		else if (pBondSmirks.getBegin() == null) {
    			begin = aam.get(rBondSmirks.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING));
    			end = aam.get(rBondSmirks.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING));
    			IBond broken = ac.getBond(begin, end);
    			ac.removeBond(broken);
    			applyAtomModification(begin, indexProduct.get(rBondSmirks.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
    			applyAtomModification(end, indexProduct.get(rBondSmirks.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
    			
    			//correct hydrogen
    			int hBegin = tools.caclulateHydrogen2(ac, begin, true);
    			if (hBegin > -1) 
    				begin.setImplicitHydrogenCount(hBegin);
//    			if (begin.getImplicitHydrogenCount() + broken.getOrder().numeric() - hBegin == 0) {
//    				begin.setImplicitHydrogenCount(hBegin);
//    			}
    			int hEnd = tools.caclulateHydrogen2(ac, end, true);
    			if (hEnd > -1) 
    				end.setImplicitHydrogenCount(hEnd);
//    			if (end.getImplicitHydrogenCount() + broken.getOrder().numeric() - hEnd == 0) {
//    				end.setImplicitHydrogenCount(hEnd);
//    			}
    		}
    		// other change s(stereo, order, atom charge,...)
    		else {
    			int orderDiff = 0;
    			begin = aam.get(pBondSmirks.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING));
    			end = aam.get(pBondSmirks.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING));
    			IBond changed = ac.getBond(begin, end);
    			changed.setAtom(begin, 0);
    			changed.setAtom(end, 1);
    			changed.setIsAromatic(pBondSmirks.isAromatic());
    			if (pBondSmirks.getOrder() != IBond.Order.UNSET) {
    				orderDiff = changed.getOrder().numeric() - pBondSmirks.getOrder().numeric();
    				changed.setOrder(pBondSmirks.getOrder());
    			}
    			changed.setStereo(pBondSmirks.getStereo());
    			applyAtomModification(begin, indexProduct.get(pBondSmirks.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
    			applyAtomModification(end, indexProduct.get(pBondSmirks.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
    		
    			//correct hydrogen
    			int hBegin = tools.caclulateHydrogen2(ac, begin, true);
    			if (hBegin > -1) 
    				begin.setImplicitHydrogenCount(hBegin);
//    			if (begin.getImplicitHydrogenCount() + orderDiff- hBegin == 0) {
//    				begin.setImplicitHydrogenCount(hBegin);
//    			}
    			int hEnd = tools.caclulateHydrogen2(ac, end, true);
    			if (hEnd > -1) 
    				end.setImplicitHydrogenCount(hEnd);
//    			if (end.getImplicitHydrogenCount() + orderDiff - hEnd == 0) {
//    				end.setImplicitHydrogenCount(hEnd);
//    			}
    		}		
    	}
    }
    
    private void applyChanges2(IAtomContainer ac, Map<IBond,IBond> changesMap, 
    		Map<Integer, List<IAtom>> aam, Map<Integer,IAtom> indexProduct) {
    	
    	Map<String, List<List<Integer>>> combinations = new HashMap<String, List<List<Integer>>>();
	
    	//si deja process ne pas ajouter ou modifier
    	for (Entry<IBond,IBond> e : changesMap.entrySet()) {
    		IBond rBondSmirks = e.getKey();
    		IBond pBondSmirks = e.getValue();
    		IAtom begin;
    		IAtom end;
    		//made
    		if (rBondSmirks.getBegin() == null) {
    			IBond newBond = new Bond();
    			int beginAam = pBondSmirks.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING);
    			int endAam = pBondSmirks.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING);
    			String id = beginAam + "-" + endAam;
    			List<IAtom> l1 = aam.get(beginAam);
    			List<IAtom> l2 = aam.get(endAam);
    			
    			if (combinations.containsKey(id)) {
    				List<List<Integer>> sols = combinations.get(id);
    				begin = l1.get(sols.get(0).get(0));
    				end = l2.get(sols.get(0).get(1));
    				sols.remove(0);
    			}
    			else {
    		    	List<List<Integer>> indexOfTheMaps = new ArrayList<List<Integer>>();
    				int index = 0;
    				List<Integer> list = new ArrayList<Integer>();
    				for (int i = 0; i < l1.size(); i++) {
    					list.add(index);
    					index++;
    				}
    				indexOfTheMaps.add(list);
    				index = 0;
    				list = new ArrayList<Integer>();
    				for (int i = 0; i < l2.size(); i++) {
    					list.add(index);
    					index++;
    				}
    				indexOfTheMaps.add(list);
    		    	List<List<Integer>> sols = new ArrayList<List<Integer>>();  
    				GenerateAllCombinationsFromMultipleLists(indexOfTheMaps, sols, 0, new ArrayList<Integer>(), false);
    				begin = l1.get(sols.get(0).get(0));
    				end = l2.get(sols.get(0).get(1));
    				sols.remove(0);
    				combinations.put(id, sols);
    			}
    			
    			if (ac.getBond(begin, end) != null) {
    				ac.removeAllElements();
    				return;
    			}
    			if (ac.getBond(begin, end) == null) {
    				newBond.setAtom(begin, 0);
        			newBond.setAtom(end, 1);
        			newBond.setIsAromatic(pBondSmirks.isAromatic());
        			newBond.setOrder(pBondSmirks.getOrder());
        			newBond.setStereo(pBondSmirks.getStereo());
        			ac.addBond(newBond);
        			applyAtomModification(begin, indexProduct.get(pBondSmirks.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
        			applyAtomModification(end, indexProduct.get(pBondSmirks.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
        			
        			//correct hydrogen
        			tools.caclulateHydrogen(ac, begin, true);
        			tools.caclulateHydrogen(ac, end, true);
        			//OLD method
        			/*
        			int hBegin = tools.caclulateHydrogen2(ac, begin, true);
        			if (hBegin > -1) 
        				begin.setImplicitHydrogenCount(hBegin);
        			int hEnd = tools.caclulateHydrogen2(ac, end, true);
        			if (hEnd > -1) 
        				end.setImplicitHydrogenCount(hEnd);
        			*/
    			}
    			
    		}
    		//broken
    		else if (pBondSmirks.getBegin() == null) {
    			int beginAam = rBondSmirks.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING);
    			int endAam = rBondSmirks.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING);
    			String id = beginAam + "-" + endAam;
    			List<IAtom> l1 = aam.get(beginAam);
    			List<IAtom> l2 = aam.get(endAam);
    			if (combinations.containsKey(id)) {
    				List<List<Integer>> sols = combinations.get(id);
    				begin = l1.get(sols.get(0).get(0));
    				end = l2.get(sols.get(0).get(1));
    				sols.remove(0);
    			}
    			else {
    		    	List<List<Integer>> indexOfTheMaps = new ArrayList<List<Integer>>();
    				int index = 0;
    				List<Integer> list = new ArrayList<Integer>();
    				for (int i = 0; i < l1.size(); i++) {
    					list.add(index);
    					index++;
    				}
    				indexOfTheMaps.add(list);
    				index = 0;
    				list = new ArrayList<Integer>();
    				for (int i = 0; i < l2.size(); i++) {
    					list.add(index);
    					index++;
    				}
    				indexOfTheMaps.add(list);
    		    	List<List<Integer>> sols = new ArrayList<List<Integer>>();  
    				GenerateAllCombinationsFromMultipleLists(indexOfTheMaps, sols, 0, new ArrayList<Integer>(), false);
    				begin = l1.get(sols.get(0).get(0));
    				end = l2.get(sols.get(0).get(1));
    				sols.remove(0);
    				combinations.put(id, sols);
    			}
    			IBond broken = ac.getBond(begin, end);
    			if (broken != null) {
    				ac.removeBond(broken);
        			applyAtomModification(begin, indexProduct.get(rBondSmirks.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
        			applyAtomModification(end, indexProduct.get(rBondSmirks.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
        			
        			//correct hydrogen
        			int hBegin = tools.caclulateHydrogen2(ac, begin, true);
        			if (hBegin > -1) 
        				begin.setImplicitHydrogenCount(hBegin);
//        			if (begin.getImplicitHydrogenCount() + broken.getOrder().numeric() - hBegin == 0) {
//        				begin.setImplicitHydrogenCount(hBegin);
//        			}
        			int hEnd = tools.caclulateHydrogen2(ac, end, true);
        			if (hEnd > -1) 
        				end.setImplicitHydrogenCount(hEnd);
//        			if (end.getImplicitHydrogenCount() + broken.getOrder().numeric() - hEnd == 0) {
//        				end.setImplicitHydrogenCount(hEnd);
//        			}
    			}
    			
    		}
    		// other change s(stereo, order, atom charge,...)
    		else {
    			int beginAam = pBondSmirks.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING);
    			int endAam = pBondSmirks.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING);
    			String id = beginAam + "-" + endAam;
    			List<IAtom> l1 = aam.get(beginAam);
    			List<IAtom> l2 = aam.get(endAam);
    			if (combinations.containsKey(id)) {
    				List<List<Integer>> sols = combinations.get(id);
    				begin = l1.get(sols.get(0).get(0));
    				end = l2.get(sols.get(0).get(1));
    				sols.remove(0);
    			}
    			else {
    		    	List<List<Integer>> indexOfTheMaps = new ArrayList<List<Integer>>();
    				int index = 0;
    				List<Integer> list = new ArrayList<Integer>();
    				for (int i = 0; i < l1.size(); i++) {
    					list.add(index);
    					index++;
    				}
    				indexOfTheMaps.add(list);
    				index = 0;
    				list = new ArrayList<Integer>();
    				for (int i = 0; i < l2.size(); i++) {
    					list.add(index);
    					index++;
    				}
    				indexOfTheMaps.add(list);
    		    	List<List<Integer>> sols = new ArrayList<List<Integer>>();  
    				GenerateAllCombinationsFromMultipleLists(indexOfTheMaps, sols, 0, new ArrayList<Integer>(), false);
    				begin = l1.get(sols.get(0).get(0));
    				end = l2.get(sols.get(0).get(1));
    				sols.remove(0);
    				combinations.put(id, sols);
    			}
    			IBond changed = ac.getBond(begin, end);
    			
    			if (changed == null) {
    				ac.removeAllElements();
    				return;
    			}
    			
    			changed.setAtom(begin, 0);
    			changed.setAtom(end, 1);
    			changed.setIsAromatic(pBondSmirks.isAromatic());
    			if (pBondSmirks.getOrder() != IBond.Order.UNSET) {
    				//orderDiff = changed.getOrder().numeric() - pBondSmirks.getOrder().numeric();
    				changed.setOrder(pBondSmirks.getOrder());
    			}
    			changed.setStereo(pBondSmirks.getStereo());
    			applyAtomModification(begin, indexProduct.get(pBondSmirks.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
    			applyAtomModification(end, indexProduct.get(pBondSmirks.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
//    		System.out.println("ccccc "+changed);
    			//correct hydrogen
    			int hBegin = tools.caclulateHydrogen2(ac, begin, true);
    			if (hBegin > -1) 
    				begin.setImplicitHydrogenCount(hBegin);
//    			if (begin.getImplicitHydrogenCount() + orderDiff- hBegin == 0) {
//    				begin.setImplicitHydrogenCount(hBegin);
//    			}
    			int hEnd = tools.caclulateHydrogen2(ac, end, true);
    			if (hEnd > -1) 
    				end.setImplicitHydrogenCount(hEnd);
//    			if (end.getImplicitHydrogenCount() + orderDiff - hEnd == 0) {
//    				end.setImplicitHydrogenCount(hEnd);
//    			}
    		}		
    	}
    }
    
    private void applyAtomModification(IAtom toModify, IAtom ref) {
    	if (ref.getCharge() != null)
    		toModify.setCharge(ref.getCharge());
    	if (ref.getExactMass() != null)
    		toModify.setExactMass(ref.getExactMass());
    	if (ref.getFormalCharge() != null)
    		toModify.setFormalCharge(ref.getFormalCharge());
    	if (ref.isAromatic() == true)
    		toModify.setIsAromatic(ref.isAromatic());
    	if ((int)ref.getProperty("bcount") != 1) {
    		if (ref.getHybridization() != null)
        		toModify.setHybridization(ref.getHybridization());
        	if (ref.getImplicitHydrogenCount() != null)
        		toModify.setImplicitHydrogenCount(ref.getImplicitHydrogenCount());
        	if (ref.getValency() != null)
        		toModify.setValency(ref.getValency());
    	}
    }


	public void setStoichiometry(int stoichiometry) {
		this.stoichiometry = stoichiometry;
	}


	public void setCheckValence(boolean checkValence) {
		this.checkValence = checkValence;
	}
    
}
