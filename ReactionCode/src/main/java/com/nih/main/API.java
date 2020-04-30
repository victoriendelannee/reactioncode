package com.nih.main;

import java.util.List;
import java.util.Map;
import java.util.Set;

import org.openscience.cdk.CDKConstants;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;


import com.nih.codes.DecodeReactionCode;
import com.nih.codes.EncodeReactionCode;
import com.nih.reaction.PseudoMolecule;
import com.nih.tools.tools;
import com.nih.transformer.Transformer;

public class API {

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
	 */
	public String encoder(IReaction reaction, boolean bondType, boolean charge, boolean hybridization, 
			boolean repetition, boolean stereochemistry) {
		String reactionCode = null;
		Map<Object, Object> properties = reaction.getProperties();
		IAtomContainerSet reactants = reaction.getReactants();
		IAtomContainerSet products = reaction.getProducts();
		try {
			for (IAtomContainer ac : reactants.atomContainers()) {
				tools.perceiveAromaticityAndflagAtomsAndBonds(ac);
				tools.attributeIDtoAtomsAndBonds(ac);
			}

			for (IAtomContainer ac : products.atomContainers()) {
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

		try {
			pm.reactionAnnotator(reaction);
			pseudoMolecule = pm.makePseudoMolecule(reactants, products);
			pseudoMolecule.setProperty(CDKConstants.TITLE, "PSEUDO-MOLECULE");
		}
		catch (Exception e) {
			System.err.println(e);
		}
		
		try {
			Map<String, Integer> numberOfRepetitions = pm.atomRepetition();
			Set<IAtom> reactionCenterAtom = pm.getReactioncenter();
			
			//make reactionCode
			//set parameters for generation of reactionCode
			EncodeReactionCode encoder = new EncodeReactionCode();
			
			//set parameters for generation of reactionCode
			encoder.setBondType(bondType);
			encoder.setCharge(charge);
			encoder.setHybridization(hybridization);
			encoder.setRepetition(repetition);
			encoder.setStereochemistry(stereochemistry);
			Map<String,String> reactionCodeMap = encoder.makeReactionCode(reactionCenterAtom, reactants, 
					pseudoMolecule, numberOfRepetitions);
			//String reactionCode = ReactionCode.reactionCodeMapToStringMultiLines(reactionCodeMap);
			reactionCode = encoder.reactionCodeMapToStringOneLine(reactionCodeMap);
		}
		catch (Exception e) {
			System.err.println(e);
		}
		return reactionCode;
	}
	
	/**
	 * Decode the reaction from the reactionCode (String) to IReacion
	 * @param reactionCode
	 * @return
	 */
	public IReaction decoder(String reactionCode) {
		DecodeReactionCode decoder = new DecodeReactionCode();
		IReaction reaction = null;
		try {
			reaction = decoder.decode(reactionCode);
			reaction.setID(0+"");
		}
		catch (Exception e){
			System.err.println(e);
		}
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
	public IReaction tranformer(String reactionCode, IAtomContainerSet reagents) throws CloneNotSupportedException, CDKException {
		IReaction transform = decoder(reactionCode);
		Transformer transformer = new Transformer();
		return transformer.applyTranform(reagents, transform);
	}
	
	/**
	 * Generate all possible products
	 * @param reactionCode
	 * @param reagents
	 * @return
	 * @throws CloneNotSupportedException
	 * @throws CDKException
	 */
	public List<IReaction> tranformer2(String reactionCode, IAtomContainerSet reagents) throws CloneNotSupportedException, CDKException {
		IReaction transform = decoder(reactionCode);
		Transformer transformer = new Transformer();
		return transformer.applyTranform2(reagents, transform);
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
	public List<IReaction> tranformer2(String reactionCode, IAtomContainerSet reagents, int stoichiometry) throws CloneNotSupportedException, CDKException {
		IReaction transform = decoder(reactionCode);
		Transformer transformer = new Transformer();
		transformer.setStoichiometry(stoichiometry);
		return transformer.applyTranform2(reagents, transform);
	}
}
