package com.nih.main;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.smiles.SmilesParser;

import com.nih.codes.DecodeReactionCode;
import com.nih.codes.EncodeReactionCode;
import com.nih.reaction.PseudoMolecule;
import com.nih.reaction.additionalConstants;
import com.nih.tools.tools;
import com.nih.transformer.Transformer;

public class API {
	Set<Integer> atomsNotInProducts;
	/**
	 * Generate Reaction code
	 * default options: bondType=true, charge=true, hybridization=false, repetition=true and stereochemistry=true
	 * alternative options: bondType=false, charge=false, hybridization=true, repetition=false and stereochemistry=false
	 * @param reaction
	 * @param bondType
	 * @param charge
	 * @param hybridization
	 * @param repetition
	 * @param stereochemistry
	 * @return
	 * @throws CDKException 
	 */
	public String encode(IReaction reaction, boolean bondType, boolean charge, boolean hybridization, 
			boolean repetition, boolean stereochemistry, boolean perceiveAromaticity) throws CDKException {
		String reactionCode = null;
		Map<Object, Object> properties = reaction.getProperties();
		IAtomContainerSet reactants = reaction.getReactants();
		IAtomContainerSet products = reaction.getProducts();
		try {
			for (IAtomContainer ac : reactants.atomContainers()) {
				if (perceiveAromaticity)
					tools.perceiveAromaticityAndflagAtomsAndBonds(ac);
				tools.attributeIDtoAtomsAndBonds(ac);
			}

			for (IAtomContainer ac : products.atomContainers()) {
				if (perceiveAromaticity)
					tools.perceiveAromaticityAndflagAtomsAndBonds(ac);
				tools.attributeIDtoAtomsAndBonds(ac);
			}
		}
		catch (Exception e) {
			System.err.println(e);
		}

		//Annotate the reaction with right constants
		PseudoMolecule pm = new PseudoMolecule();
		IAtomContainer pseudoMolecule = null;

			pm.reactionAnnotator(reaction);
			
			pseudoMolecule = pm.makePseudoMolecule(reactants, products);
			pseudoMolecule.setProperty(CDKConstants.TITLE, "PSEUDO-MOLECULE");

		//TODO update main.java and uncomment try
		try {
			Map<String, Integer> numberOfRepetitions = pm.atomRepetition();
			Set<IAtom> reactionCenterAtom = pm.getReactioncenter();
		//make reactionCode
		//set parameters for generation of reactionCode
			EncodeReactionCode encoder = new EncodeReactionCode(new HashSet<Integer>());
			
			//set parameters for generation of reactionCode
			encoder.setBondType(bondType);
			encoder.setCharge(charge);
			encoder.setHybridization(hybridization);
			encoder.setRepetition(repetition);
			encoder.setStereochemistry(stereochemistry);
			Map<String,String> reactionCodeMap = encoder.generateReactionCode(reactionCenterAtom, reactants, 
					pseudoMolecule, numberOfRepetitions);
			//String reactionCode = ReactionCode.reactionCodeMapToStringMultiLines(reactionCodeMap);
			reactionCode = encoder.reactionCodeMapToStringOneLine(reactionCodeMap);
			atomsNotInProducts = encoder.getAtomInReactionCenterNotInProducts();
		}
		catch (Exception e) {
			System.err.println(e);
		}
		return reactionCode;
	}	
	
	/**
	 * Generate Reaction code and restore the balance of unbalanced reactions
	 * default options: bondType=true, charge=true, hybridization=false, repetition=true and stereochemistry=true
	 * alternative options: bondType=false, charge=false, hybridization=true, repetition=false and stereochemistry=false
	 * @param reaction
	 * @param bondType
	 * @param charge
	 * @param hybridization
	 * @param repetition
	 * @param stereochemistry
	 * @return
	 */
	public String encodeAndRebalance(IReaction reaction, boolean bondType, boolean charge, boolean hybridization, 
			boolean repetition, boolean stereochemistry, boolean perceiveAromaticity) {

		String reactionCode = null;
		try {
			reactionCode = encode(reaction, bondType, charge, hybridization,repetition, stereochemistry, false);
		} catch (CDKException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		if (!reactionCode.contains("A:"))
			return reactionCode;
		IReaction decoded;
		try {
			decoded = decode(reactionCode, false, true,  atomsNotInProducts);
			removeProperties(decoded.getReactants());
			removeProperties(decoded.getProducts());
			return encode(decoded, bondType, charge, hybridization, 
					 repetition, stereochemistry, false);
		} catch (CDKException | CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
		
		
	}
	
	private void removeProperties(IAtomContainerSet set) {
		for (IAtomContainer ac : set.atomContainers()) {
			removeProperties(ac);
		}
	}
	
	private void removeProperties(IAtomContainer ac) {
		for (IAtom a : ac.atoms()) {
			int aam = a.getProperty(CDKConstants.ATOM_ATOM_MAPPING);
			a.setProperties(new HashMap());
			a.setProperty(CDKConstants.ATOM_ATOM_MAPPING, aam);
		}
		for (IBond b : ac.bonds()) {
			b.setProperties(new HashMap());
		}
	}
	
	/**
	 * Decode the reaction from the reactionCode (String) to IReaction
	 * @param reactionCode
	 * @param kekulize
	 * @param correctProducts
	 * @param notInProducts
	 * @return
	 * @throws CDKException
	 * @throws CloneNotSupportedException
	 */
	private IReaction decode(String reactionCode, boolean kekulize, boolean correctProducts, Set<Integer> notInProducts) throws CDKException, CloneNotSupportedException {
		DecodeReactionCode decoder = new DecodeReactionCode(notInProducts);
		decoder.kekulize(kekulize);
		decoder.setCorrectProducts(correctProducts);
		IReaction reaction = null; 
		try {
			reaction = decoder.decode(reactionCode);
			reaction.setID(0+"");
		}
		catch (Exception e){
			System.err.println(e);
		}
		reaction.setProperty("errors", decoder.getErrors());
		return reaction;
	}
	
	/**
	 * Decode the reaction from the reactionCode (String) to IReaction
	 * @param reactionCode
	 * @param kekulize
	 * @param correctProducts
	 * @return
	 * @throws CDKException
	 * @throws CloneNotSupportedException
	 */
	public IReaction decode(String reactionCode, boolean kekulize, boolean correctProducts) throws CDKException, CloneNotSupportedException {
		DecodeReactionCode decoder = new DecodeReactionCode();
		decoder.kekulize(kekulize);
		decoder.setCorrectProducts(correctProducts);
		IReaction reaction = null; 
		try {
			reaction = decoder.decode(reactionCode);
			reaction.setID(0+"");
		}
		catch (Exception e){
			System.err.println(e);
		}
		reaction.setProperty("errors", decoder.getErrors());
		return reaction;
	}
	
	/**
	 * Decode the reaction from the reactionCode (String) to IReaction
	 * @param reactionCode
	 * @return
	 * @throws CloneNotSupportedException 
	 * @throws CDKException 
	 */
	public IReaction decode(String reactionCode) throws CDKException, CloneNotSupportedException {
		DecodeReactionCode decoder = new DecodeReactionCode();
		IReaction reaction = null;
		try {
			reaction = decoder.decode(reactionCode);
			reaction.setID(0+"");
		}
		catch (Exception e){
			System.err.println(e);
		}
		reaction.setProperty("errors", decoder.getErrors());
		return reaction;
	}
	
	/**
	 * Generate One solution (best match)
	 * @param reactionCode
	 * @param reagents
	 * @return
	 * @throws CloneNotSupportedException
	 * @throws CDKException
	 */
	public IReaction tranform(String reactionCode, IAtomContainerSet reagents) throws CloneNotSupportedException, CDKException {
		IReaction transform = decode(reactionCode);
		Transformer transformer = new Transformer();
		transformer.setStoichiometry(1);
		return transformer.transform(reagents, transform).get(0);
	}
	
	/**
	 * Generate all possible products
	 * @param reactionCode
	 * @param reagents
	 * @return
	 * @throws CloneNotSupportedException
	 * @throws CDKException
	 */
	public List<IReaction> tranform2(String reactionCode, IAtomContainerSet reagents) throws CloneNotSupportedException, CDKException {
		IReaction transform = decode(reactionCode);
		Transformer transformer = new Transformer();
		return transformer.transform(reagents, transform);
	}
	
	/**
	 * Generate all possible products
	 * @param reactionCode
	 * @param reagents
	 * @param stoichiometry (number of times a same reactant can react on different site)
	 * @return
	 * @throws CloneNotSupportedException
	 * @throws CDKException
	 */
	public List<IReaction> tranform2(String reactionCode, IAtomContainerSet reagents, int stoichiometry) throws CloneNotSupportedException, CDKException {
		IReaction transform = decode(reactionCode);
		Transformer transformer = new Transformer();
		transformer.setStoichiometry(stoichiometry);
		return transformer.transform(reagents, transform);
	}
	
	/**
	 * @param reaction
	 * @param atomAtomMapping
	 * @return
	 * @throws CDKException 
	 */
	public String makeSMILES(IReaction reaction, boolean atomAtomMapping) throws CDKException {
		return tools.makeSmiles(reaction, atomAtomMapping); 
	}
	
	/**
	 * @param set
	 * @param atomAtomMapping
	 * @return
	 * @throws CDKException
	 */
	public String makeSMILES(IAtomContainerSet set, boolean atomAtomMapping) throws CDKException {
		String smi = "";
		for (IAtomContainer ac : set.atomContainers()) {
			smi +=makeSMILES(ac, atomAtomMapping) + ".";
		}
		smi = smi.substring(0, smi.length() - 1);
		return smi;
	}
	
	/**
	 * @param ac
	 * @param atomAtomMapping
	 * @return
	 * @throws CDKException
	 */
	public String makeSMILES(IAtomContainer ac, boolean atomAtomMapping) throws CDKException {
		return tools.makeSmiles(ac, true, atomAtomMapping, false) ;
	}
}
