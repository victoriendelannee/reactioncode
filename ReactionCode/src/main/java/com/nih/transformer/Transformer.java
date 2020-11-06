package com.nih.transformer;

import static com.nih.reaction.additionalConstants.BOND_CHANGE_INFORMATION;
import static com.nih.reaction.additionalConstants.BOND_MADE;
import static com.nih.reaction.additionalConstants.BOND_CLEAVED;
import static com.nih.reaction.additionalConstants.BOND_ORDER_CHANGE;
import static com.nih.reaction.additionalConstants.BOND_STEREO_CHANGE;
import static org.openscience.cdk.silent.SilentChemObjectBuilder.getInstance;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.openscience.cdk.silent.Atom;
import org.openscience.cdk.silent.Bond;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.fingerprint.CircularFingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.isomorphism.matchers.Expr;
import org.openscience.cdk.isomorphism.matchers.Expr.Type;
import org.openscience.cdk.isomorphism.matchers.QueryAtom;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.stereo.DoubleBondStereochemistry;
import org.openscience.cdk.stereo.TetrahedralChirality;

import com.nih.fragments.FragmentUtils;
import com.nih.tools.tools;

public class Transformer {

	int stoichiometry = 1;
	boolean checkValence = true;
	boolean allPatternsHaveToMatch = true;
	List<int[]> bondStereoChange = new ArrayList<int[]>();
	List<BitSet> productFP = new ArrayList<BitSet>();
	// store molecule, which doen't match with at least one pattern
	List<IAtomContainer> reagentsList = new ArrayList<IAtomContainer>();
	//will change if at least one bond is UNSET
	boolean kekulized = true;

	/**
	 * @param reactants
	 * @param reactionPattern
	 * @return
	 */
	public List<IReaction> transform(IAtomContainerSet reactants, IReaction reactionPattern) {
		IAtomContainerSet products = getInstance().newInstance(IAtomContainerSet.class);
		int cpt = 1;
		for (IAtomContainer ac : reactants.atomContainers()) {
			try {
				Kekulization.kekulize(ac);
			} catch (CDKException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			for (IAtom atom : ac.atoms()) {
				atom.removeProperty(CDKConstants.ATOM_ATOM_MAPPING);
				atom.setID(cpt + "");
				// atom.setProperty(CDKConstants.ATOM_ATOM_MAPPING, cpt);
				cpt++;
			}
			products.addAtomContainer(ac);
		}
		
		return transform(products, reactionPattern, reactants);
	}

	/**
	 * @param reactantsSmiles
	 * @param reactionPattern
	 * @return
	 * @throws InvalidSmilesException
	 */
	public List<IReaction> transform(String reactantsSmiles, IReaction reactionPattern) throws InvalidSmilesException {
		SmilesParser sp = new SmilesParser(getInstance());
		sp.kekulise(true);
		String[] reactantsSmi = reactantsSmiles.split("\\.");

		IAtomContainerSet reactants = getInstance().newInstance(IAtomContainerSet.class);
		IAtomContainerSet products = getInstance().newInstance(IAtomContainerSet.class);

		// prepare products container to apply the reaction and and index bonds
		int cpt = 1;
		for (String reactant : reactantsSmi) {
			IAtomContainer rea = sp.parseSmiles(reactant);
			for (int i = 0; i < rea.getAtomCount(); i++) {
				IAtom atom = rea.getAtom(i);
				atom.setID(cpt + "");
				// atom.setProperty(CDKConstants.ATOM_ATOM_MAPPING, cpt);
				cpt++;
			}
			reactants.addAtomContainer(rea);
			products.addAtomContainer(rea);
		}
		return transform(products, reactionPattern, reactants);
	}

	/**
	 * @param oProducts
	 * @param reactionPattern
	 * @param reactants
	 * @return
	 */
	private List<IReaction> transform(IAtomContainerSet oProducts, IReaction reactionPattern,
			IAtomContainerSet reactants) {

		// save AAM savemapping[index] = AAM integer
		// int[] savemapping = saveMapping(aggregateProducts);
		// reset list
		productFP = new ArrayList<BitSet>();

		List<IReaction> results = new ArrayList<IReaction>();
		// get all possible reactants combinations matching with smarts
		List<List<Mappings>> ptrnMolAssociation = findPatternMoleculeAssociation(reactants,
				reactionPattern.getReactants());
		// if it returns null it means that all patterns haven't matched (it's related
		// to the global variable allPatternsHaveToMatch)
		if (ptrnMolAssociation == null)
			return results;
		// remove reagents from reactant
		IAtomContainerSet reagents = getReagents(reactants, oProducts);

		IAtomContainer aggregateProducts = aggregateIAtomContainerSet(oProducts);
		// reset Visited Flag (used by CDK when using SMARTS pAttern)
		resetVisitedFlag(aggregateProducts);

		// index of all possible mappings for all reactants
		List<List<Integer>> indexOfTheMaps = new ArrayList<List<Integer>>();
		List<List<Map<IAtom, IAtom>>> combinationList = new ArrayList<List<Map<IAtom, IAtom>>>();
		// Each Mapping has a unique ID id order to differentiate them. This map is used
		// to extract the right mapping in combination list
		// where each sublist contains the Mapping objects related to one Pattern (ex
		// one SMARTS).
		// We'll need to combine at least one mapping from each sublist in order to
		// match with all Pattern define in the ReactionCOde or SMIRKS
		Map<Integer, Integer> indexesOfMappingInIndexOfTheMaps = new HashMap<Integer, Integer>();
		int mappingID = 0;
		for (List<Mappings> combination : ptrnMolAssociation) {
			List<Map<IAtom, IAtom>> result = new ArrayList<Map<IAtom, IAtom>>();
			for (Mappings mapping : combination) {
				// convert Iterable to List
				// involved in the reaction center)
				mapping.uniqueAtoms().toAtomMap().forEach(result::add);
			}
			// get the index of all possible mapping in the list
			List<Integer> list = new ArrayList<Integer>();
			for (int i = 0; i < result.size(); i++) {
				list.add(mappingID);
				indexesOfMappingInIndexOfTheMaps.put(mappingID, i);
				mappingID++;
			}
			indexOfTheMaps.add(list);
			combinationList.add(result);
		}

		// get all possible combination
		List<List<Integer>> sols = new ArrayList<List<Integer>>();
		GenerateAllCombinationsFromMultipleLists(indexOfTheMaps, sols, 0, new ArrayList<Integer>(), false);

		List<IAtomContainer> reactantToAdd = new ArrayList<IAtomContainer>();
		boolean validMappingSolution = false;
		for (int i = 0; i < sols.size(); i++) {
			// reset mapping
			resetMappingAndBondChangeInfo(aggregateProducts);
			List<Integer> solution = sols.get(i);
			// attribute atom atom mapping to atom in aggregatesproducts
			// (which is a IAtomContainer containing all reactants and will be modified to
			// generate the products)
			for (int j = 0; j < solution.size(); j++) {
				mappingID = solution.get(j);
				int ind = indexesOfMappingInIndexOfTheMaps.get(mappingID);
				if (j > 0) {
					// same mapping ID, so another reactant has to be added (stoichiometry > 1)
					if (mappingID == solution.get(j - 1)) {
						reactantToAdd.add(
								ptrnMolAssociation.get(j).get(ind).getQuery().getProperty("correspondingReactant"));
					}
				}
				validMappingSolution = attributeAtomAtomMapping(aggregateProducts, combinationList.get(j).get(ind));
				if (!validMappingSolution) break;
			}
			
			if (!validMappingSolution) continue;


			IAtomContainer aggregateProductsCopy = null;
			try {
				aggregateProductsCopy = aggregateProducts.clone();
			} catch (CloneNotSupportedException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}

			for (IAtomContainer ac : reactantToAdd) {
				try {
					aggregateProductsCopy.add(ac.clone());
				} catch (CloneNotSupportedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}

			try {
				bondStereoChange = new ArrayList<int[]>();

				IAtomContainer newAggregateProducts = applyTransform(reactionPattern, aggregateProductsCopy,
						aggregateProducts);

				IAtomContainerSet reactantsCopy = SilentChemObjectBuilder.getInstance()
						.newInstance(IAtomContainerSet.class);
				for (IAtomContainer ac : reactants.atomContainers()) {
					reactantsCopy.addAtomContainer(ac.clone());
				}
				annotateStereochange(reactantsCopy, newAggregateProducts);

				IAtomContainerSet products = FragmentUtils.makeAtomContainerSet(newAggregateProducts);

				if (checkValence) {
					if (!tools.checkProductValidity(products)) {
						continue;
					}
				}
				
				if (!checkProductDuplicate(products)) {
					System.err.println("duplicate");
					continue;
				}

				IReaction reaction = SilentChemObjectBuilder.getInstance().newInstance(IReaction.class);
				reaction.setReactants(reactantsCopy);
				for (IAtomContainer ac : reagents.atomContainers()) {
					reaction.addAgent(ac);
				}
				reaction.setProducts(products);
				results.add(reaction);
				// convert explicit hydrogens to implicit
				convertExplicitToImplicitHydrogens(reaction);

				// correct reagents
				finalizeReaction(reaction);

			} catch (CloneNotSupportedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return results;
	}

	private void finalizeReaction(IReaction reaction) {
		List<IAtomContainer> toRemove = new ArrayList<IAtomContainer>();

		IAtomContainerSet products = reaction.getProducts();
		IAtomContainerSet reactants = reaction.getReactants();

		for (IAtomContainer ac : products.atomContainers()) {
			boolean hasReactionCenter = false;
			for (IAtom atom : ac.atoms()) {
				if (atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING) != null) {
					hasReactionCenter = true;
				}
			}
			if (!hasReactionCenter)
				toRemove.add(ac);
		}
		for (IAtomContainer ac : toRemove) {
			products.removeAtomContainer(ac);
		}

		toRemove.clear();

		for (IAtomContainer ac : reactants.atomContainers()) {
			boolean hasReactionCenter = false;
			for (IAtom atom : ac.atoms()) {
				if (atom.getProperty(BOND_CHANGE_INFORMATION) != null) {
					hasReactionCenter = true;
				}
			}
			if (!hasReactionCenter)
				toRemove.add(ac);
		}
		for (IAtomContainer ac : toRemove) {
			reaction.addAgent(ac);
			reactants.removeAtomContainer(ac);
		}

		// reattribuate mapping
		mappingReattribution(reaction);
	}

	/**
	 * @param set
	 * @return
	 */
	private IAtomContainer aggregateIAtomContainerSet(IAtomContainerSet set) {
		IAtomContainer aggregate = SilentChemObjectBuilder.getInstance().newAtomContainer();

		for (IAtomContainer ac : set.atomContainers()) {
			aggregate.add(ac);
		}
		return aggregate;
	}

	/**
	 * Return the Reagents (Molecules which are the same in reactants and products)
	 * 
	 * @param reactants
	 * @return
	 */
	private IAtomContainerSet getReagents(IAtomContainerSet reactants, IAtomContainerSet products) {
		IAtomContainerSet reagents = SilentChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		for (IAtomContainer reagent : reagentsList) {
			reactants.removeAtomContainer(reagent);
			products.removeAtomContainer(reagent);
			reagents.addAtomContainer(reagent);
		}
		return reagents;
	}

	/**
	 * Annotate the stereo bond changes
	 * 
	 * @param reactants
	 * @param aggregateProduct
	 */
	private void annotateStereochange(IAtomContainerSet reactants, IAtomContainer aggregateProduct) {
		IAtomContainer aggregateReactants = SilentChemObjectBuilder.getInstance().newAtomContainer();
		for (IAtomContainer ac : reactants.atomContainers()) {
			aggregateReactants.add(ac);
		}
		for (int[] arr : bondStereoChange) {
			IAtom beginR = getAtomByAam(aggregateReactants, arr[0]);
			IAtom endR = getAtomByAam(aggregateReactants, arr[1]);
			IBond bondR = aggregateReactants.getBond(beginR, endR);

			IAtom beginP = getAtomByAam(aggregateProduct, arr[0]);
			IAtom endP = getAtomByAam(aggregateProduct, arr[1]);
			IBond bondP = aggregateProduct.getBond(beginP, endP);

			if (!bondR.getStereo().equals(bondP.getStereo())) {
				bondR.setProperty(BOND_CHANGE_INFORMATION, BOND_STEREO_CHANGE);
				bondR.getBegin().setProperty(BOND_CHANGE_INFORMATION, BOND_STEREO_CHANGE);
				bondR.getEnd().setProperty(BOND_CHANGE_INFORMATION, BOND_STEREO_CHANGE);
				bondP.setProperty(BOND_CHANGE_INFORMATION, BOND_STEREO_CHANGE);
				bondP.getBegin().setProperty(BOND_CHANGE_INFORMATION, BOND_STEREO_CHANGE);
				bondP.getEnd().setProperty(BOND_CHANGE_INFORMATION, BOND_STEREO_CHANGE);
			}
		}
	}

	/**
	 * Check if the product has been already generated
	 * 
	 * @param aggregateProducts
	 * @return
	 */
	private boolean checkProductDuplicate(IAtomContainerSet set) {
		if (!set.isEmpty()) {
			// Calculate FP and check if the product was already generated

			for (IAtomContainer ac : set.atomContainers()) {
				try {
					CircularFingerprinter fp = new CircularFingerprinter();
					fp.calculate(ac);
					IBitFingerprint bfp = fp.getBitFingerprint(ac);
					if (productFP.contains(bfp.asBitSet()))
						return false;
				} catch (CDKException e) {
					// TODO Auto-generated catch block
					// e.printStackTrace();
					return true;
				}
			}
			
		} else
			return false;
		return true;
	}

	private void resetVisitedFlag(IAtomContainer ac) {
		for (IAtom atom : ac.atoms()) {
			atom.setFlag(CDKConstants.VISITED, false);
		}
	}

	private int[] saveMapping(IAtomContainer ac) {
		int[] save = new int[ac.getAtomCount()];
		for (int i = 0; i < ac.getAtomCount(); i++) {
			IAtom atom = ac.getAtom(i);
			save[i] = atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING);
		}
		return save;
	}

	private void resetMappingAndBondChangeInfo(IAtomContainer ac) {
		for (IAtom atom : ac.atoms()) {
			atom.removeProperty(CDKConstants.ATOM_ATOM_MAPPING);
			atom.removeProperty(BOND_CHANGE_INFORMATION);
		}
		for (IBond bond : ac.bonds()) {
			bond.removeProperty(BOND_CHANGE_INFORMATION);
			if (bond.getOrder().equals(IBond.Order.UNSET))
				kekulized = false;
		}
	}

	private boolean testMapping(IAtomContainerSet set) {
		for (IAtomContainer ac : set.atomContainers()) {
			if (!testMapping(ac))
				return false;
		}
		return true;
	}

	private boolean testMapping(IAtomContainer ac) {
		for (IAtom atom : ac.atoms()) {
			if (atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING) != null)
				return true;
		}
		return false;
	}

	private Mappings getMapping(IAtomContainer reactantPattern, IAtomContainer ac) {
		Pattern ptrn = Pattern.findSubstructure(reactantPattern);
		Mappings m = ptrn.matchAll(ac);
		if (m.count() == 0)
			return null;
		else
			return m;
	}

	private List<List<Mappings>> findPatternMoleculeAssociation(IAtomContainerSet reactants,
			IAtomContainerSet reactantPatterns) {
		List<List<Mappings>> result = new ArrayList<List<Mappings>>();

		Mappings[][] possibilities = new Mappings[reactants.getAtomContainerCount()][reactantPatterns
				.getAtomContainerCount()];
		for (int i = 0; i < reactants.getAtomContainerCount(); i++) {
			IAtomContainer ac = reactants.getAtomContainer(i);
			boolean matched = false;
			// match with each SMARTS pattern in reactants side
			for (int j = 0; j < reactantPatterns.getAtomContainerCount(); j++) {
				IAtomContainer reactantPattern = reactantPatterns.getAtomContainer(j);
				Mappings mappings = getMapping(reactantPattern, ac);
				if (mappings == null) {
					possibilities[i][j] = null;
				} else {
					possibilities[i][j] = mappings;
					matched = true;
				}
			}
			if (!matched)
				reagentsList.add(ac);
		}
		// get pattern association ex pattern 1 matches with r1 r2 and ptrn 2 with r1 r2
		// r3 result=[[r1,r2],[r1,r2,r3]]
		for (int i = 0; i < reactantPatterns.getAtomContainerCount(); i++) {
			List<Mappings> mapping = new ArrayList<Mappings>();
			for (int j = 0; j < reactants.getAtomContainerCount(); j++) {
				Mappings m = possibilities[j][i];
				if (m != null) {
					m.getQuery().setProperty("correspondingReactant", reactants.getAtomContainer(j));
					mapping.add(m);
				}
			}
			if (mapping.isEmpty() && allPatternsHaveToMatch) {
				return null;
			}
			result.add(mapping);
		}
		return result;
	}

	/**
	 * for each reactant, get the mappings of all matching SMARTS store the
	 * individual mapping results (null if no mapping) SMARTS1 SMARTS2 SMARTS3 ac1
	 * L1 null L2 ac2 null L7 L3 ac3 L4 L5 L6
	 * 
	 * @param reactants
	 * @param smirks
	 * @return
	 */
	private Mappings[][] matchAllPatternsWithAllReactants(IAtomContainerSet reactants,
			IAtomContainerSet reactantPatterns) {
		Mappings[][] possibilities = new Mappings[reactants.getAtomContainerCount()][reactantPatterns
				.getAtomContainerCount()];
		for (int i = 0; i < reactants.getAtomContainerCount(); i++) {
			IAtomContainer ac = reactants.getAtomContainer(i);

			boolean matched = false;
			// match with each SMARTS pattern in reactants side
			for (int j = 0; j < reactantPatterns.getAtomContainerCount(); j++) {
				IAtomContainer reactantPattern = reactantPatterns.getAtomContainer(j);

				Mappings mappings = getMapping(reactantPattern, ac);
				if (mappings == null) {
					possibilities[i][j] = null;
				} else {
					possibilities[i][j] = mappings;
					matched = true;
				}
			}
			if (!matched)
				reagentsList.add(ac);
		}
		// remove molecules which don't match with any patterns
		Mappings[][] possibilities2 = new Mappings[reactants.getAtomContainerCount()
				- reagentsList.size()][reactantPatterns.getAtomContainerCount()];
		int index = 0;
		for (int i = 0; i < reactants.getAtomContainerCount(); i++) {
			IAtomContainer ac = reactants.getAtomContainer(i);
			if (!reagentsList.contains(ac)) {
				for (int j = 0; j < reactantPatterns.getAtomContainerCount(); j++) {
					possibilities2[index][j] = possibilities[i][j];
				}
				index++;
			}
		}
		return possibilities2;
	}

	/**
	 * Return all possible combinations (those combining all SMARTS, ex L1,L5,L3;
	 * L1,L7,L6; L4,L7,L2) (cf tab in findCompatibleSmarts) A valid combination is a
	 * combination, where each reactant is associated with one SMARTS result =
	 * [[L1,L5,L3,] [L1,L7,L6], [L4,L7,L2]]
	 * 
	 * @param possibilities
	 * @return
	 */
	private List<List<Mappings>> getCompatibleMappingAssociation(Mappings[][] possibilities) {
		int numberOfPatterns = possibilities[0].length;
		List<List<Mappings>> combinations = new ArrayList<List<Mappings>>();
		for (int i = 0; i < numberOfPatterns; i++) {
			if (i == 0) {
				List<Mappings> combination = new ArrayList<Mappings>();
				combinations.add(combination);
			}
			for (int j = 0; j < possibilities.length; j++) {
				if (possibilities[j][i] != null) {
					for (List<Mappings> combination : combinations) {
						combination.add(possibilities[j][i]);
					}
				}
			}
		}
		return combinations;
	}

	// take hydrogen count(explicit and Implicit), charge bond order, bond
	// made/broken, aromaticity (bond + atom), stereo
	private IAtomContainer applyTransform(IReaction reactionPattern, IAtomContainer aggregateProducts,
			IAtomContainer aggregateReactants) {
		// aggregate reactants and product in smirks
		QueryAtomContainer reactantsPattern = new QueryAtomContainer(SilentChemObjectBuilder.getInstance());
		QueryAtomContainer productsPattern = new QueryAtomContainer(SilentChemObjectBuilder.getInstance());
		for (IAtomContainer ac : reactionPattern.getReactants().atomContainers()) {
			reactantsPattern.add(ac);
		}
		for (IAtomContainer ac : reactionPattern.getProducts().atomContainers()) {
			productsPattern.add(ac);
		}

		// register reactants mapping
		List<Integer> reactantRecords = registerMappings(reactantsPattern);
		List<Integer> productRecords = registerMappings(productsPattern);

		// check empty bond
		List<IBond> bondsToDel = new ArrayList<IBond>();
		for (IBond bond : aggregateProducts.bonds()) {
			if (bond.getBegin() == null || bond.getEnd() == null)
				bondsToDel.add(bond);
		}
		for (IBond bond : bondsToDel) {
			aggregateProducts.removeBond(bond);
		}

		for (IBond bond : productsPattern.bonds()) {
			QueryBond qb = (QueryBond) bond;
			boolean addBond = false;
			QueryAtom begin = (QueryAtom) bond.getBegin();
			QueryAtom end = (QueryAtom) bond.getEnd();
			IAtom configBegin = null;
			IAtom configEnd = null;
			// all atoms without mapping are considered absent in reactants and are created
			if (begin.getProperty("done") == null) {
				if (begin.getProperty(CDKConstants.ATOM_ATOM_MAPPING) == null) {
					// test if multiple possible atom type
					IAtom atom = addNewAtom(begin, productsPattern, aggregateProducts);
					configBegin = atom;
					addBond = true;
					// flag
					begin.setProperty("done", true);
					begin.setProperty("configAtom", atom);
				} else {
					// get atom in aggregate product and configure it
					if (reactantRecords.contains(begin.getProperty(CDKConstants.ATOM_ATOM_MAPPING))) {
						int aam = begin.getProperty(CDKConstants.ATOM_ATOM_MAPPING);
						configBegin = configAtom((QueryAtom) reactantsPattern.getAtom(reactantRecords.indexOf(aam)),
								begin, getAtomByAam(aggregateProducts, aam), aggregateProducts);
						begin.setProperty("configAtom", configBegin);
					} else {
						configBegin = configAtom(aggregateProducts, begin);
						begin.setProperty("configAtom", configBegin);
					}
				}
			} else
				configBegin = begin.getProperty("configAtom");

			if (end.getProperty("done") == null) {
				if (end.getProperty(CDKConstants.ATOM_ATOM_MAPPING) == null) {
					// test if multiple possible atom type
					IAtom atom = addNewAtom(end, productsPattern, aggregateProducts);
					configEnd = atom;
					addBond = true;
					// flag
					end.setProperty("done", true);
					end.setProperty("configAtom", atom);
				} else {
					// get atom in aggregate product and configure it
					if (reactantRecords.contains(end.getProperty(CDKConstants.ATOM_ATOM_MAPPING))) {
						int aam = end.getProperty(CDKConstants.ATOM_ATOM_MAPPING);
						configEnd = configAtom((QueryAtom) reactantsPattern.getAtom(reactantRecords.indexOf(aam)), end,
								getAtomByAam(aggregateProducts, aam), aggregateProducts);
						end.setProperty("configAtom", configEnd);
					} else {
						configEnd = configAtom(aggregateProducts, end);
						end.setProperty("configAtom", configEnd);
					}
				}
			} else
				configEnd = end.getProperty("configAtom");
			if (addBond == true) {
				IBond formedBond = makeBond(aggregateProducts, configBegin, configEnd, qb);
				formedBond.setProperty(BOND_CHANGE_INFORMATION, BOND_MADE);
				formedBond.getBegin().setProperty(BOND_CHANGE_INFORMATION, BOND_MADE);
				formedBond.getEnd().setProperty(BOND_CHANGE_INFORMATION, BOND_MADE);
			} else {
				if (begin.getProperty(CDKConstants.ATOM_ATOM_MAPPING) == null
						|| end.getProperty(CDKConstants.ATOM_ATOM_MAPPING) == null) {
					IBond formedBond = makeBond(aggregateProducts, configBegin, configEnd, qb);
					formedBond.setProperty(BOND_CHANGE_INFORMATION, BOND_MADE);
					formedBond.getBegin().setProperty(BOND_CHANGE_INFORMATION, BOND_MADE);
					formedBond.getEnd().setProperty(BOND_CHANGE_INFORMATION, BOND_MADE);
				}
				// modify bond if necessary
				else if (reactantsPattern.getBond(
						reactantsPattern
								.getAtom(reactantRecords.indexOf(begin.getProperty(CDKConstants.ATOM_ATOM_MAPPING))),
						reactantsPattern.getAtom(
								reactantRecords.indexOf(end.getProperty(CDKConstants.ATOM_ATOM_MAPPING)))) != null) {
					QueryBond qReactant = (QueryBond) reactantsPattern.getBond(
							reactantsPattern.getAtom(
									reactantRecords.indexOf(begin.getProperty(CDKConstants.ATOM_ATOM_MAPPING))),
							reactantsPattern
									.getAtom(reactantRecords.indexOf(end.getProperty(CDKConstants.ATOM_ATOM_MAPPING))));
					IBond aggegateBond = aggregateProducts.getBond(
							getAtomByAam(aggregateProducts, begin.getProperty(CDKConstants.ATOM_ATOM_MAPPING)),
							getAtomByAam(aggregateProducts, end.getProperty(CDKConstants.ATOM_ATOM_MAPPING)));

					configBond(qReactant, qb, aggegateBond);
					aggegateBond.setProperty(BOND_CHANGE_INFORMATION, BOND_ORDER_CHANGE);
					aggegateBond.getBegin().setProperty(BOND_CHANGE_INFORMATION, BOND_ORDER_CHANGE);
					aggegateBond.getEnd().setProperty(BOND_CHANGE_INFORMATION, BOND_ORDER_CHANGE);
				}
				// bond made
				else {
					IBond formedBond = makeBond(aggregateProducts, configBegin, configEnd, qb);
					formedBond.setProperty(BOND_CHANGE_INFORMATION, BOND_MADE);
					formedBond.getBegin().setProperty(BOND_CHANGE_INFORMATION, BOND_MADE);
					formedBond.getEnd().setProperty(BOND_CHANGE_INFORMATION, BOND_MADE);
				}
			}

		}

		for (IBond rpbond : reactantsPattern.bonds()) {
			IAtom rpBegin = rpbond.getBegin();
			IAtom rpEnd = rpbond.getEnd();

			IBond ppBond = productsPattern.getBond(
					productsPattern.getAtom(productRecords.indexOf(rpBegin.getProperty(CDKConstants.ATOM_ATOM_MAPPING))),
					productsPattern.getAtom(productRecords.indexOf(rpEnd.getProperty(CDKConstants.ATOM_ATOM_MAPPING))));
			IBond pBond = aggregateProducts.getBond(
					getAtomByAam(aggregateProducts, rpBegin.getProperty(CDKConstants.ATOM_ATOM_MAPPING)),
					getAtomByAam(aggregateProducts, rpEnd.getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
			IBond rBond = aggregateReactants.getBond(
					getAtomByAam(aggregateReactants, rpBegin.getProperty(CDKConstants.ATOM_ATOM_MAPPING)),
					getAtomByAam(aggregateReactants, rpEnd.getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
			// means cleaved
			if (ppBond == null) {
				rBond.setProperty(BOND_CHANGE_INFORMATION, BOND_CLEAVED);
				rBond.getBegin().setProperty(BOND_CHANGE_INFORMATION, BOND_CLEAVED);
				rBond.getEnd().setProperty(BOND_CHANGE_INFORMATION, BOND_CLEAVED);
				pBond.getBegin().setProperty(BOND_CHANGE_INFORMATION, BOND_CLEAVED);
				pBond.getEnd().setProperty(BOND_CHANGE_INFORMATION, BOND_CLEAVED);
				aggregateProducts.removeBond(pBond);
			} else {
				if (pBond.getProperty(BOND_CHANGE_INFORMATION) != null) {
					rBond.setProperty(BOND_CHANGE_INFORMATION, pBond.getProperty(BOND_CHANGE_INFORMATION));
					rBond.getBegin().setProperty(BOND_CHANGE_INFORMATION, pBond.getProperty(BOND_CHANGE_INFORMATION));
					rBond.getEnd().setProperty(BOND_CHANGE_INFORMATION, pBond.getProperty(BOND_CHANGE_INFORMATION));
				}
			}

		}
		for (IAtom rpatom : reactantsPattern.atoms()) {
			if (rpatom.getProperty(BOND_CHANGE_INFORMATION) == null) {
				IAtom pAtom = getAtomByAam(aggregateProducts, rpatom.getProperty(CDKConstants.ATOM_ATOM_MAPPING));
				IAtom rAtom = getAtomByAam(aggregateReactants, rpatom.getProperty(CDKConstants.ATOM_ATOM_MAPPING));
				if (pAtom != null) {
					rAtom.setProperty(BOND_CHANGE_INFORMATION, pAtom.getProperty(BOND_CHANGE_INFORMATION));
				}
				else {
					rAtom.setProperty(BOND_CHANGE_INFORMATION, BOND_CLEAVED);
				}
			}

		}
		// finalize stereo
		finalize(aggregateProducts, productsPattern);

		return aggregateProducts;
	}

	private IAtom addNewAtom(QueryAtom qAtom, QueryAtomContainer qac, IAtomContainer products) {
		Expr.Type[] testUniqueElement = { Expr.Type.ALIPHATIC_ELEMENT, Expr.Type.AROMATIC_ELEMENT, Expr.Type.ELEMENT };
		Map<Type, List<Integer>> exprs = qAtom.getExpression().toMap();
		boolean matches = false;
		for (Expr.Type type : testUniqueElement) {
			if (exprs.containsKey(type)) {
				if (matches) {
					System.err.println("Can not manage yet more than one element type");
					System.exit(0);
				}
				if (exprs.get(type).size() > 1) {
					System.err.println("Can not manage yet more than one element type");
					System.exit(0);
				} else {
					// calculate delta by getting min in list
					IAtom atom = new Atom(exprs.get(type).get(0));
					products.addAtom(atom);
					if (exprs.containsKey(Expr.Type.FORMAL_CHARGE)) {
						if (exprs.get(Expr.Type.FORMAL_CHARGE).size() > 1) {
							System.err.println("Can not manage yet more than one possible charge ");
							System.exit(0);
						} else
							atom.setFormalCharge(exprs.get(Expr.Type.FORMAL_CHARGE).get(0));
					}
					if (exprs.containsKey(Expr.Type.IMPL_H_COUNT)) {
						atom.setImplicitHydrogenCount(exprs.get(Expr.Type.IMPL_H_COUNT).get(0));
					}
					// get explicit hydrogen and deduce implicit number
					if (exprs.containsKey(Expr.Type.TOTAL_H_COUNT)) {
						int explicitHydrogen = exprs.get(Expr.Type.TOTAL_H_COUNT).get(0);
						int hCount = 0;
						for (IAtom atom2 : qac.getConnectedAtomsList(atom)) {
							QueryAtom qa = (QueryAtom) atom2;
							Map<Type, List<Integer>> exprs2 = qa.getExpression().toMap();
							for (Expr.Type type2 : testUniqueElement) {
								if (exprs2.containsKey(type)) {
									if (matches) {
										System.err.println("Can not manage yet more than one element type");
										System.exit(0);
									}
									if (exprs.get(type2).size() > 1) {
										System.err.println("Can not manage yet more than one element type");
										System.exit(0);
									} else {
										if (exprs.get(type2).get(0) == 1)
											hCount++;
									}
								}
							}
						}
						// create missing explicit hydrogen
						if (hCount != explicitHydrogen) {
							for (int i = 0; i < explicitHydrogen - hCount; i++) {
								Atom newHydrogen = new Atom("H");
								products.addAtom(newHydrogen);
								Bond newBond = new Bond();
								newBond.setAtom(atom, 0);
								newBond.setAtom(newHydrogen, 1);
								newBond.setOrder(IBond.Order.SINGLE);
								products.addBond(newBond);
							}
						}
						// atom.setImplicitHydrogenCount(exprs.get(Expr.Type.TOTAL_H_COUNT).get(0) -
						// hCount);
					}
					if (exprs.containsKey(Expr.Type.IS_AROMATIC) || type.equals(Expr.Type.AROMATIC_ELEMENT))
						atom.setIsAromatic(true);
					else
						atom.setIsAromatic(false);
					if (exprs.containsKey(Expr.Type.IS_IN_RING) || exprs.containsKey(Expr.Type.RING_COUNT)
							|| exprs.containsKey(Expr.Type.RING_SIZE))
						atom.setIsInRing(true);
					else
						atom.setIsInRing(false);
					if (exprs.containsKey(Expr.Type.VALENCE)) {
						if (exprs.get(Expr.Type.VALENCE).size() > 1) {
							System.err.println("Can not manage yet more than one possible valence ");
							System.exit(0);
						} else
							atom.setFormalCharge(exprs.get(Expr.Type.VALENCE).get(0));
					}
					return atom;
				}
			}
		}
		return null;
	}

	private List<Integer> registerMappings(IAtomContainer ac) {
		List<Integer> records = new ArrayList<Integer>();
		for (int i = 0; i < ac.getAtomCount(); i++) {
			IAtom atom = ac.getAtom(i);
			if (atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING) != null)
				records.add(atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING));
			else
				records.add(null);
			// reset custom properties
			atom.removeProperty("configAtom");
			atom.removeProperty("done");
		}
		return records;
	}

	private void configBond(QueryBond source, QueryBond target, IBond toModify) {
		Map<Type, List<Integer>> exprs1 = source.getExpression().toMap();
		Map<Type, List<Integer>> exprs2 = target.getExpression().toMap();
		for (Entry<Type, List<Integer>> entry : exprs2.entrySet()) {
			Expr.Type type = entry.getKey();
			List<Integer> val2 = entry.getValue();
			if (exprs1.containsKey(type)) {
				List<Integer> val1 = exprs1.get(type);
				if (type.equals(Expr.Type.ALIPHATIC_ORDER)) {
					int delta = findMin(val2) - findMin(val1);
					int order = toModify.getOrder().numeric() + delta;
					switch (order) {
					case 1:
						toModify.setOrder(IBond.Order.SINGLE);
						break;
					case 2:
						toModify.setOrder(IBond.Order.DOUBLE);
						break;
					case 3:
						toModify.setOrder(IBond.Order.TRIPLE);
						break;
					case 4:
						toModify.setOrder(IBond.Order.QUADRUPLE);
						break;
					}

				}
			}
			if (type.equals(Expr.Type.IS_AROMATIC) || type.equals(Expr.Type.AROMATIC_ELEMENT)) {
				toModify.setIsAromatic(true);
				toModify.getBegin().setIsAromatic(true);
				toModify.getEnd().setIsAromatic(true);
				if (toModify.getOrder().equals(IBond.Order.UNSET))
					toModify.setOrder(IBond.Order.SINGLE);
			} else
				toModify.setIsAromatic(false);
			if (type.equals(Expr.Type.IS_IN_RING) || type.equals(Expr.Type.RING_COUNT)
					|| type.equals(Expr.Type.RING_SIZE))
				toModify.setIsInRing(true);
			else
				toModify.setIsInRing(false);
			if (type.equals(Expr.Type.STEREOCHEMISTRY)) {
				int stereo = val2.get(0);
				switch (stereo) {
				case 0x2f:
					toModify.setStereo(IBond.Stereo.UP);
					break;
				case 0x5c:
					toModify.setStereo(IBond.Stereo.DOWN);
					break;
				case 0x9a:
					toModify.setStereo(IBond.Stereo.UP_INVERTED);
					break;
				case 0x9b:
					toModify.setStereo(IBond.Stereo.DOWN_INVERTED);
					break;
				case 0x9c:
					toModify.setStereo(IBond.Stereo.E);
					break;
				case 0x9d:
					toModify.setStereo(IBond.Stereo.Z);
					break;
				case 0x9e:
					toModify.setStereo(IBond.Stereo.E_OR_Z);
					break;
				case 0x9f:
					toModify.setStereo(IBond.Stereo.UP_OR_DOWN);
					break;
				case 0x8f:
					toModify.setStereo(IBond.Stereo.UP_OR_DOWN_INVERTED);
					break;
				case 0x00:
					toModify.setStereo(IBond.Stereo.NONE);
					break;
				}
			}
		}
	}

	private IBond makeBond(IAtomContainer ac, IAtom begin, IAtom end, QueryBond ref) {
		IBond newBond = new Bond();
		newBond.setAtom(begin, 0);
		newBond.setAtom(end, 1);
		newBond.setIsAromatic(false);
		newBond.setOrder(IBond.Order.SINGLE);
		Map<Type, List<Integer>> exprs = ref.getExpression().toMap();
		if (exprs.containsKey(Expr.Type.SINGLE_OR_AROMATIC)) {
			newBond.setOrder(IBond.Order.SINGLE);
		} else if (exprs.containsKey(Expr.Type.DOUBLE_OR_AROMATIC)) {
			newBond.setOrder(IBond.Order.DOUBLE);
		}
		if (exprs.containsKey(Expr.Type.ALIPHATIC_ORDER)) {
			int order = exprs.get(Expr.Type.ALIPHATIC_ORDER).get(0);
			switch (order) {
			case 1:
				newBond.setOrder(IBond.Order.SINGLE);
				break;
			case 2:
				newBond.setOrder(IBond.Order.DOUBLE);
				break;
			case 3:
				newBond.setOrder(IBond.Order.TRIPLE);
				break;
			case 4:
				newBond.setOrder(IBond.Order.QUADRUPLE);
				break;
			}

		}
		if (exprs.containsKey(Expr.Type.IS_AROMATIC)) {
			newBond.setIsAromatic(true);
			newBond.getBegin().setIsAromatic(true);
			newBond.getEnd().setIsAromatic(true);
			newBond.setOrder(IBond.Order.SINGLE);
		}
		if (exprs.containsKey(Expr.Type.STEREOCHEMISTRY)) {
			int stereo = exprs.get(Expr.Type.STEREOCHEMISTRY).get(0);
			switch (stereo) {
			case 0x2f:
				newBond.setStereo(IBond.Stereo.UP);
				break;
			case 0x5c:
				newBond.setStereo(IBond.Stereo.DOWN);
				break;
			}
		}
		if (exprs.containsKey(Expr.Type.IS_IN_RING) || exprs.containsKey(Expr.Type.RING_COUNT)
				|| exprs.containsKey(Expr.Type.RING_SIZE))
			newBond.setIsInRing(true);
		else
			newBond.setIsInRing(false);
		ac.addBond(newBond);
		return newBond;
	}

	private IAtom configAtom(QueryAtom source, QueryAtom target, IAtom toModify, IAtomContainer ac) {
		Map<Type, List<Integer>> exprs1 = source.getExpression().toMap();
		Map<Type, List<Integer>> exprs2 = target.getExpression().toMap();
		for (Entry<Type, List<Integer>> entry : exprs2.entrySet()) {
			Expr.Type type = entry.getKey();
			List<Integer> val2 = entry.getValue();
			if (exprs1.containsKey(type)) {
				List<Integer> val1 = exprs1.get(type);
				if (type.equals(Expr.Type.FORMAL_CHARGE)) {
					int delta = findMin(val2) - findMin(val1);
					int charge = toModify.getFormalCharge() + delta;
					toModify.setFormalCharge(charge);
				} else if (type.equals(Expr.Type.IMPL_H_COUNT)) {
					int delta = findMin(val2) - findMin(val1);
					int hCount = toModify.getImplicitHydrogenCount() + delta;
					toModify.setImplicitHydrogenCount(hCount);
				} else if (type.equals(Expr.Type.TOTAL_H_COUNT)) {
					int delta = findMin(val2) - findMin(val1);
					int hCount = toModify.getImplicitHydrogenCount();
					// add hydrogen
					if (delta > 0) {
						toModify.setImplicitHydrogenCount(hCount);
					}
					// remove hydrogen
					else if (delta < 0) {
						int cpt = 0;
						// first remove implicit hydrogen and then explicit if necessary
						while (cpt != delta) {
							if (hCount + delta >= 0) {
								hCount--;
								toModify.setImplicitHydrogenCount(hCount);
							} else {
								for (IBond bond : ac.getConnectedBondsList(toModify)) {
									IAtom other = bond.getOther(toModify);
									if (other.getSymbol().equals("H")) {
										ac.removeBond(bond);
										ac.removeAtom(other);
									}
								}
							}
							cpt--;
						}
					}
				} else if (type.equals(Expr.Type.VALENCE)) {
					int delta = findMin(val2) - findMin(val1);
					int valence = toModify.getValency() + delta;
					toModify.setValency(valence);

				}
			} else {
				if (type.equals(Expr.Type.FORMAL_CHARGE)) {
					if (val2.size() > 1) {
						System.err.println("Can not manage yet more than one FORMAL_CHARGE type");
						System.exit(0);
					} else
						toModify.setFormalCharge(val2.get(0));
				} else if (type.equals(Expr.Type.IMPL_H_COUNT)) {
					if (val2.size() > 1) {
						System.err.println("Can not manage yet more than one IMPL_H_COUNT type");
						System.exit(0);
					} else
						toModify.setImplicitHydrogenCount(val2.get(0));
				} else if (type.equals(Expr.Type.TOTAL_H_COUNT)) {
					if (val2.size() > 1) {
						System.err.println("Can not manage yet more than one TOTAL_H_COUNT type");
						System.exit(0);
					} else
						toModify.setImplicitHydrogenCount(val2.get(0));
				} else if (type.equals(Expr.Type.VALENCE)) {
					if (val2.size() > 1) {
						System.err.println("Can not manage yet more than one VALENCE type");
						System.exit(0);
					} else
						toModify.setValency(val2.get(0));

				}
			}
			if (type.equals(Expr.Type.IS_AROMATIC) || type.equals(Expr.Type.AROMATIC_ELEMENT))
				toModify.setIsAromatic(true);
			else
				toModify.setIsAromatic(false);
			if (type.equals(Expr.Type.IS_IN_RING) || type.equals(Expr.Type.RING_COUNT)
					|| type.equals(Expr.Type.RING_SIZE))
				toModify.setIsInRing(true);
			else
				toModify.setIsInRing(false);
		}
		return toModify;
	}

	private IAtom configAtom(IAtomContainer ac, QueryAtom ref) {
		IAtom target = getAtomByAam(ac, ref.getProperty(CDKConstants.ATOM_ATOM_MAPPING));
		Expr.Type[] testUniqueElement = { Expr.Type.ALIPHATIC_ELEMENT, Expr.Type.AROMATIC_ELEMENT, Expr.Type.ELEMENT };
		Map<Type, List<Integer>> exprs = ref.getExpression().toMap();
		boolean matches = false;
		for (Expr.Type type : testUniqueElement) {
			if (exprs.containsKey(type)) {
				if (matches) {
					System.err.println("Can not manage yet more than one element type");
					System.exit(0);
				}
				if (exprs.get(type).size() > 1) {
					System.err.println("Can not manage yet more than one element type");
					System.exit(0);
				} else {
					// calculate delta by getting min in list
					if (exprs.containsKey(Expr.Type.FORMAL_CHARGE)) {
						if (exprs.get(Expr.Type.FORMAL_CHARGE).size() > 1) {
							System.err.println("Can not manage yet more than one possible charge ");
							System.exit(0);
						} else
							target.setFormalCharge(exprs.get(Expr.Type.FORMAL_CHARGE).get(0));
					}
					if (exprs.containsKey(Expr.Type.IMPL_H_COUNT))
						target.setImplicitHydrogenCount(exprs.get(Expr.Type.IMPL_H_COUNT).get(0));
					// get explicit hydrogen and deduce implicit number
					if (exprs.containsKey(Expr.Type.TOTAL_H_COUNT)) {
						int hCount = 0;
						for (IAtom atom2 : ac.getConnectedAtomsList(target)) {
							QueryAtom qa = (QueryAtom) atom2;
							Map<Type, List<Integer>> exprs2 = qa.getExpression().toMap();
							for (Expr.Type type2 : testUniqueElement) {
								if (exprs.containsKey(type)) {
									if (matches) {
										System.err.println("Can not manage yet more than one element type");
										System.exit(0);
									}
									if (exprs.get(type2).size() > 1) {
										System.err.println("Can not manage yet more than one element type");
										System.exit(0);
									} else {
										if (exprs.get(type).get(0) == 1)
											hCount++;
									}
								}
							}
							target.setImplicitHydrogenCount(exprs.get(Expr.Type.IMPL_H_COUNT).get(0) - hCount);
						}
						if (exprs.containsKey(Expr.Type.IS_AROMATIC) || type.equals(Expr.Type.AROMATIC_ELEMENT))
							target.setIsAromatic(true);
						else
							target.setIsAromatic(false);
						if (exprs.containsKey(Expr.Type.IS_IN_RING) || exprs.containsKey(Expr.Type.RING_COUNT)
								|| exprs.containsKey(Expr.Type.RING_SIZE))
							target.setIsInRing(true);
						else
							target.setIsInRing(false);
						if (exprs.containsKey(Expr.Type.VALENCE)) {
							if (exprs.get(Expr.Type.VALENCE).size() > 1) {
								System.err.println("Can not manage yet more than one possible valence ");
								System.exit(0);
							} else
								target.setFormalCharge(exprs.get(Expr.Type.VALENCE).get(0));
						}
					}
					return target;
				}
			}
		}
		return target;
	}

	private void finalize(IAtomContainer aggregate, QueryAtomContainer productsSmirks) {
		Map<IChemObject, IStereoElement> savedStereo = new HashMap<IChemObject, IStereoElement>();
		List<IStereoElement> stereo = new ArrayList<IStereoElement>();
		for (IStereoElement se : aggregate.stereoElements()) {
			savedStereo.put(se.getFocus(), se);
			stereo.add(se);
		}
		aggregate.setStereoElements(new ArrayList<IStereoElement>());
		for (IStereoElement se : productsSmirks.stereoElements()) {
			if (se instanceof TetrahedralChirality) {
				TetrahedralChirality th = (TetrahedralChirality) se;
				IAtom chiral = getAtomByAam(aggregate, th.getChiralAtom().getProperty(CDKConstants.ATOM_ATOM_MAPPING));
				IAtom[] ligands = new IAtom[4];
				for (int i = 0; i < th.getLigands().length; i++) {
					IAtom atom = th.getLigands()[i];
					ligands[i] = getAtomByAam(aggregate, atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING));
				}
				if (savedStereo.containsKey(chiral))
					stereo.remove(savedStereo.get(chiral));
				stereo.add(new TetrahedralChirality(chiral, ligands, th.getConfig()));
				// comment because it destroyed stereo outside of the pattern
				// aggregate.addStereoElement(new TetrahedralChirality(chiral, ligands,
				// th.getConfig()));

				// get the bond with the stereo change
				for (IBond bond : aggregate.getConnectedBondsList(chiral)) {
					if (bond.getStereo() != null) {
						if (bond.getStereo() != IBond.Stereo.NONE) {
							int[] arr = new int[2];
							arr[0] = bond.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING);
							arr[1] = bond.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING);
							bondStereoChange.add(arr);
						}
					}
				}
			} else if (se instanceof DoubleBondStereochemistry) {
				DoubleBondStereochemistry db = (DoubleBondStereochemistry) se;
				IBond focus = aggregate.getBond(
						getAtomByAam(aggregate, db.getFocus().getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING)),
						getAtomByAam(aggregate, db.getFocus().getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
				IBond[] ligands = new IBond[2];
				for (int i = 0; i < db.getBonds().length; i++) {
					IBond bond = db.getBonds()[i];
					ligands[i] = aggregate.getBond(
							getAtomByAam(aggregate, bond.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING)),
							getAtomByAam(aggregate, bond.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
				}
				if (savedStereo.containsKey(focus))
					stereo.remove(savedStereo.get(focus));
				stereo.add(new DoubleBondStereochemistry(focus, ligands, db.getConfig()));
				// comment because it destroyed stereo outside of the pattern
				// aggregate.addStereoElement(new DoubleBondStereochemistry(focus, ligands,
				// db.getConfig()));

				// get the bond with the stereo change
				int[] arr = new int[2];
				arr[0] = focus.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING);
				arr[1] = focus.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING);
				bondStereoChange.add(arr);
			}
		}
		aggregate.setStereoElements(stereo);

		// correct hydrogen for atom involve in a bond made or broken
		for (IAtom atom : aggregate.atoms()) {
			if (atom.getProperty(BOND_CHANGE_INFORMATION) != null) {
				//System.out.println(atom.getProperty(BOND_CHANGE_INFORMATION) + " " +atom); 
				if (atom.getProperty(BOND_CHANGE_INFORMATION).equals(BOND_MADE)
						|| atom.getProperty(BOND_CHANGE_INFORMATION).equals(BOND_CLEAVED)
						|| atom.getProperty(BOND_CHANGE_INFORMATION).equals(BOND_ORDER_CHANGE)) {
	    			tools.caclulateHydrogen(aggregate, atom, kekulized);
	    			//OLD Method
	    			/*
					int hCount = tools.caclulateHydrogen(aggregate, atom, kekulized);
					if (hCount > -1)
						atom.setImplicitHydrogenCount(hCount);
					*/
				}
			}
		}
	}

	private IAtom getAtomByID(IAtomContainerSet set, String id) {
		for (IAtomContainer ac : set.atomContainers()) {
			IAtom atom = getAtomByID(ac, id);
			if (atom != null)
				return atom;
		}
		return null;
	}

	private IAtom getAtomByID(IAtomContainer ac, String id) {
		for (IAtom atom : ac.atoms()) {
			if (atom.getID() != null) {
				if (atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING).equals(id))
					return atom;
			}
		}
		return null;
	}

	private IAtom getAtomByAam(IAtomContainerSet set, int aam) {
		for (IAtomContainer ac : set.atomContainers()) {
			IAtom atom = getAtomByAam(ac, aam);
			if (atom != null)
				return atom;
		}
		return null;
	}

	private IAtom getAtomByAam(IAtomContainer ac, int aam) {
		for (IAtom atom : ac.atoms()) {
			if (atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING) != null) {
				if (atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING).equals(aam))
					return atom;
			}
		}
		return null;
	}

	public void mappingReattribution(IReaction reaction) {
		for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
			for (IAtom atom : ac.atoms()) {
				atom.setProperty(CDKConstants.ATOM_ATOM_MAPPING, Integer.parseInt(atom.getID()));
			}
		}
		for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
			for (IAtom atom : ac.atoms()) {
				atom.setProperty(CDKConstants.ATOM_ATOM_MAPPING, Integer.parseInt(atom.getID()));
			}
		}
	}
	// the new method take the id and assign it as AAM
	/*
	 * public void mappingReattribution(IReaction reaction) { //old mapping new
	 * mapping Map<Integer,Integer> aam = new TreeMap<Integer,Integer>();
	 * 
	 * for (IAtomContainer ac : reaction.getReactants().atomContainers()) { for
	 * (IAtom atom : ac.atoms()) {
	 * aam.put(atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING), -1); } } int
	 * counter = 1; for (Entry<Integer,Integer> e : aam.entrySet()) { int i =
	 * e.getKey(); aam.put(i, counter); counter++; } for (IAtomContainer ac :
	 * reaction.getReactants().atomContainers()) { for (IAtom atom : ac.atoms()) {
	 * int newAam = aam.get(atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING));
	 * atom.setProperty(CDKConstants.ATOM_ATOM_MAPPING, newAam); } } for
	 * (IAtomContainer ac : reaction.getProducts().atomContainers()) { for (IAtom
	 * atom : ac.atoms()) { Integer oldAam =
	 * atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING); if (oldAam == null) {
	 * atom.setProperty(CDKConstants.ATOM_ATOM_MAPPING, counter); counter++; } else
	 * if (!aam.containsKey(oldAam)){
	 * atom.setProperty(CDKConstants.ATOM_ATOM_MAPPING, counter); counter++; } else
	 * atom.setProperty(CDKConstants.ATOM_ATOM_MAPPING, aam.get(oldAam)); } } }
	 */

	// function to find minimum value in an unsorted
	// list in Java using Collection
	public Integer findMin(List<Integer> list) {

		// check list is empty or not
		if (list == null || list.size() == 0) {
			return Integer.MAX_VALUE;
		}

		// create a new list to avoid modification
		// in the original list
		List<Integer> sortedlist = new ArrayList<>(list);

		// sort list in natural order
		Collections.sort(sortedlist);

		// first element in the sorted list
		// would be minimum
		return sortedlist.get(0);
	}

	// function return maximum value in an unsorted
	// list in Java using Collection
	public Integer findMax(List<Integer> list) {

		// check list is empty or not
		if (list == null || list.size() == 0) {
			return Integer.MIN_VALUE;
		}

		// create a new list to avoid modification
		// in the original list
		List<Integer> sortedlist = new ArrayList<>(list);

		// sort list in natural order
		Collections.sort(sortedlist);

		// last element in the sorted list would be maximum
		return sortedlist.get(sortedlist.size() - 1);
	}

	private boolean attributeAtomAtomMapping(IAtomContainer ac, Map<IAtom, IAtom> aam) {
		// key:query atom in SMARTS value: atom in AC
		for (Entry<IAtom, IAtom> e : aam.entrySet()) {
			IAtom queryAtom = e.getKey();
			IAtom targetAtom = e.getValue();
			//the mapping algorithm matched with the same atom, when it should match with different atom
			//so one atom in reactant is never match (O.OC>>OC.O both oxygen should match 2 different 
			//atoms but they can match the same atom with CDK )
			if (targetAtom.getProperty(CDKConstants.ATOM_ATOM_MAPPING) != null)
				return false;
			
			if (queryAtom.getProperty(CDKConstants.ATOM_ATOM_MAPPING) == null)
				targetAtom.setProperty("noMapping", true);
			else
				targetAtom.setProperty(CDKConstants.ATOM_ATOM_MAPPING,
						(int) queryAtom.getProperty(CDKConstants.ATOM_ATOM_MAPPING));
		}
		return true;
	}

	/**
	 * Generate combination from multiple lists X: [1, 2] Y: [3, 4, 5, 6] Result:
	 * [[1, 3], [1, 4], [1, 5], [1, 6], [2, 3], [2, 4], [2, 5], [2, 6]] if parameter
	 * intraCombinations is true, generate combination between the lists and for
	 * each combination generate all permutations inside the permuted list X: [1, 2]
	 * [3, 4, 5, 6] Result: [[1, 3], [1, 4], [1, 5], [1, 6], [2, 3], [2, 4], [2, 5],
	 * [2, 6], [1, 4, 3], [1, 5, 3], [1, 5, 4], [1, 6, 3], [1, 6, 4], [1, 6, 5], [2,
	 * 4, 3], [2, 5, 3], [2, 5, 4], [2, 6, 3], [2, 6, 4], [2, 6, 5], [2, 1, 3], [2,
	 * 1, 4], [2, 1, 5], [2, 1, 6], [1, 5, 3, 4], [1, 6, 3, 4], [1, 6, 3, 5], [1, 6,
	 * 4, 5], [2, 5, 3, 4], [2, 6, 3, 4], [2, 6, 3, 5], [2, 6, 4, 5], [2, 1, 4, 3],
	 * [2, 1, 5, 3], [2, 1, 5, 4], [2, 1, 6, 3], [2, 1, 6, 4], [2, 1, 6, 5], [1, 6,
	 * 3, 4, 5], [2, 6, 3, 4, 5], [2, 1, 5, 3, 4], [2, 1, 6, 3, 4], [2, 1, 6, 3, 5],
	 * [2, 1, 6, 4, 5], [2, 1, 6, 3, 4, 5]]
	 * https://stackoverflow.com/questions/17192796/generate-all-combinations-from-multiple-lists
	 * https://javahungry.blogspot.com/2014/02/algorithm-for-combinations-of-string-java-code-with-example.html
	 * https://www.geeksforgeeks.org/print-all-possible-combinations-of-r-elements-in-a-given-array-of-size-n/
	 * 
	 * @param Lists
	 * @param result
	 * @param depth
	 * @param current
	 */
	private void GenerateAllCombinationsFromMultipleLists(List<List<Integer>> Lists, List<List<Integer>> result,
			int depth, List<Integer> current, boolean intraCombinations) {
		if (depth == Lists.size()) {
			result.add(new ArrayList<Integer>(current));
			return;
		}

		for (int i = 0; i < Lists.get(depth).size(); i++) {
			int t = Lists.get(depth).get(i);

			List<Integer> temp = new ArrayList<Integer>(current);
			temp.add(t);

			GenerateAllCombinationsFromMultipleLists(Lists, result, depth + 1, temp, intraCombinations);

			// generate all combination inside a list and return a list which contains lists
			// of different combination
			// ex sublist [3, 4, 5] -> combinationResultsOfTheSublist [[3, 4, 5], [3, 4],
			// [3, 5], [3], [4, 5], [4], [5]]
			// limit to 5 times (stoichiometry)
			if (intraCombinations == true && i > 0 && i < stoichiometry) {
				List<Integer> sublist = Lists.get(depth).subList(0, i);
				Integer[] sublistArray = sublist.toArray(new Integer[sublist.size()]);
				List<List<Integer>> combinationResultsOfTheSublist = new ArrayList<List<Integer>>();
				generateAllCombinationsInsideAList(sublistArray, combinationResultsOfTheSublist);

				// generate all combinations of the previously generated results with the other
				// list in Lists
				for (int j = 0; j < combinationResultsOfTheSublist.size(); j++) {
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
	private void recursive_combinations(List<Integer> combination, int ndx, Integer[] elems,
			List<List<Integer>> result) {
		if (ndx == elems.length) {

			// (reached end of list after selecting/not selecting)
			if (!combination.isEmpty())
				result.add(new ArrayList<Integer>(combination));

		} else {

			// (include element at ndx)
			combination.add(elems[ndx]);
			recursive_combinations(combination, ndx + 1, elems, result);

			// (don't include element at ndx)
			combination.remove(elems[ndx]);
			recursive_combinations(combination, ndx + 1, elems, result);

		}
	}

	/**
	 * @param reaction
	 * @param aam
	 * @return
	 * @throws CDKException
	 */
	public String makeSmiles(IReaction reaction, boolean aam) throws CDKException {
		SmilesGenerator sg;
		if (aam == false)
			sg = new SmilesGenerator(SmiFlavor.Stereo | SmiFlavor.UseAromaticSymbols);
		else
			sg = new SmilesGenerator(SmiFlavor.Stereo | SmiFlavor.AtomAtomMap | SmiFlavor.UseAromaticSymbols);
		return sg.create(reaction);
	}

	/**
	 * @param reaction
	 * @param aam
	 * @return
	 * @throws CDKException
	 */
	public String makeSmiles(IAtomContainer ac, boolean aam) throws CDKException {
		SmilesGenerator sg;
		if (aam == false)
			sg = new SmilesGenerator(SmiFlavor.Stereo | SmiFlavor.UseAromaticSymbols);
		else
			sg = new SmilesGenerator(SmiFlavor.Stereo | SmiFlavor.AtomAtomMap | SmiFlavor.UseAromaticSymbols);
		return sg.create(ac);
	}

	public void convertExplicitToImplicitHydrogens(IReaction reaction) {
		for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
			convertExplicitToImplicitHydrogens(ac);
		}
		for (IAtomContainer ac : reaction.getAgents().atomContainers()) {
			convertExplicitToImplicitHydrogens(ac);
		}
		for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
			convertExplicitToImplicitHydrogens(ac);
		}
	}

	/**
	 * Returns The number of explicit hydrogens for a given IAtom.
	 * 
	 * @param atomContainer
	 * @param atom
	 * @return The number of explicit hydrogens on the given IAtom.
	 */
	public void convertExplicitToImplicitHydrogens(IAtomContainer atomContainer) {
		List<IAtom> toDelete = new ArrayList<IAtom>();
		for (IAtom atom : atomContainer.atoms()) {
			if (atom.getSymbol().equals("H"))
				continue;
			int hCount = atom.getImplicitHydrogenCount();
			for (IAtom connectedAtom : atomContainer.getConnectedAtomsList(atom)) {
				if (connectedAtom.getSymbol().equals("H")) {
					for (IBond bond : atomContainer.getConnectedBondsList(connectedAtom)) {
						if (bond.getStereo() == null) {
							atomContainer.removeBond(bond);
							toDelete.add(connectedAtom);
							hCount++;
						} else {
							if (bond.getStereo() == IBond.Stereo.NONE) {
								atomContainer.removeBond(bond);
								toDelete.add(connectedAtom);
								hCount++;
							}
						}
					}
				}
			}
			atom.setImplicitHydrogenCount(hCount);
		}
		for (IAtom atom : toDelete) {
			atomContainer.removeAtom(atom);
		}
	}

	public int getStoichiometry() {
		return stoichiometry;
	}

	public void setStoichiometry(int stoichiometry) {
		this.stoichiometry = stoichiometry;
	}

	public boolean isCheckValence() {
		return checkValence;
	}

	public void setCheckValence(boolean checkValence) {
		this.checkValence = checkValence;
	}

}
