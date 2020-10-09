package com.nih.tools;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IReaction;
//import org.openscience.cdk.interfaces.IRingSet;
//import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;

import com.nih.reaction.additionalConstants;

public class tools {

	/**
	 * @param ac
	 * @throws CDKException
	 */
	public static void attributeIDtoAtomsAndBonds(IAtomContainer ac) throws CDKException {
		// attribute ID
		for (IBond b : ac.bonds()) {
			List<String> ids = new ArrayList<String>();
			ids.add(b.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING).toString());
			ids.add(b.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING).toString());
			Collections.sort(ids);
			b.setID(ids.get(0) + "-" + ids.get(1));
			b.setProperty("aromatic", b.isAromatic());
		}
		for (IAtom a : ac.atoms()) {
			a.setID(a.getProperty(CDKConstants.ATOM_ATOM_MAPPING).toString());
		}
	}

	/**
	 * @param ac
	 * @throws CDKException
	 */
	public static void perceiveAromaticityAndflagAtomsAndBonds(IAtomContainer ac) throws CDKException {
		ElectronDonation model = ElectronDonation.daylight();
		CycleFinder cycles = Cycles.all();
		Aromaticity aromaticity = new Aromaticity(model, cycles);

		// perceive aromaticity
		aromaticity.apply(ac);

		/*
		 * //flag bond and atom in ring AllRingsFinder arf = new AllRingsFinder();
		 * IRingSet r = arf.findAllRings(ac); // manually flag ring bonds for
		 * (IAtomContainer ar : r.atomContainers()) { for (IBond bond : ar.bonds())
		 * bond.setFlag(CDKConstants.ISINRING, true); }
		 */
	}


	/**
	 * @param ac
	 * @param maxDepth
	 * @param kekulized
	 */
	public static void addMissingHydrogen(IAtomContainer ac, int maxDepth, boolean kekulized) {
		for (IAtom atom : ac.atoms()) {
			if (maxDepth > -1) {
				if ((int) atom.getProperty("depth") == maxDepth)
					continue;
			}
			if (atom.getSymbol().equals("R")) {
				atom.setImplicitHydrogenCount(0);
				continue;
			}
			//double bondOrderSum = ac.getBondOrderSum(atom);
			double bondOrderSum = 0;
			boolean isAromatic = false;
			// System.out.println(atom.getSymbol() + " "+ac.getBondOrderSum(atom) + " "
			// +ac.getConnectedBondsCount(atom));
			for (IBond bond : ac.getConnectedBondsList(atom)) {
				IBond.Order order = bond.getOrder();
				// System.out.println(order + " " + order.numeric() + " " +bond.isAromatic());
				if (order.equals(IBond.Order.UNSET) && bond.isAromatic()){
					bondOrderSum += 1;
					isAromatic = true;
				}
				else if (!order.equals(IBond.Order.UNSET) && bond.isAromatic()){
					if (kekulized) {
						bondOrderSum += order.numeric();
					}
					else {
						bondOrderSum += 1;
						isAromatic = true;
					}
				}
				else if (!order.equals(IBond.Order.UNSET) && !bond.isAromatic())
					bondOrderSum += order.numeric();
			}
			//adjust aromatic bond score temp fix if c(c)c score should be 4 if cc 3 and 3 for cn(C)n (will be adjusted below)
			if (isAromatic) bondOrderSum += 1;

			int hcount = atom.getImplicitHydrogenCount() == CDKConstants.UNSET ? 0 : atom.getImplicitHydrogenCount();
			int charge = atom.getFormalCharge() == CDKConstants.UNSET ? 0 : atom.getFormalCharge();

			double valence;
			// calculate the valence
			// TODO integrate min and max valence
//			if (atom.getBondOrderSum() != null) {
//				valence = atom.getBondOrderSum();
//			}
//			else {
			valence = ElementCalculation.calculateValence(atom.getSymbol());
//			}

			if (bondOrderSum - charge + hcount != valence) {
				int diff = (int) (valence - (bondOrderSum - charge + hcount));
				if (diff < 0 && isAromatic) diff += 1;
				int implicitHydrogen2 = atom.getImplicitHydrogenCount() + diff;
				int implicitHydrogen = implicitHydrogen2 > 0 ? implicitHydrogen2 : 0;
				atom.setImplicitHydrogenCount(implicitHydrogen);
			}
		}
	}

	public static boolean checkProductValidity(IAtomContainerSet acset) {
		// check if generated product is valid (only check modified atoms)
		for (IAtomContainer ac : acset.atomContainers()) {
			boolean isValencyValid = isValencyOfMoleculeValid(ac);
			if (!isValencyValid) {
				return false;
			}

			// used previously to check all the atoms of all molecules
			/*
			 * boolean isValencyValid = tools.isValencyOfMoleculeValid(ac); if
			 * (!isValencyValid){ valid = false; }
			 */
		}
		return true;
	}
	
	/**
	 * @param ac
	 * @return
	 */
	public static boolean isValencyOfMoleculeValid(IAtomContainer ac) {
		boolean mustReturn = false;
		try {
			Kekulization.kekulize(ac);
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			// e.printStackTrace();
			mustReturn = true;
			// return true;
		}
		
		/*SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Canonical); 
		 try {
		 System.out.println("smi " +sg.create(ac)); 
		 } catch (CDKException e) { 
			 // TODOAuto-generated catch block
			 e.printStackTrace();
			 }*/
		 

		if (mustReturn)
			return true;
		for (IAtom atom : ac.atoms()) {
			//test only atom in reaction center
			if (atom.getProperty(additionalConstants.BOND_CHANGE_INFORMATION) == null) 
				continue;
			if (atom.getSymbol().equals("R")) {
				atom.setImplicitHydrogenCount(0);
				continue;
			}
			if (ac.getConnectedBondsCount(atom) == 0) {
				continue;
			}
			double bondOrderSum = ac.getBondOrderSum(atom);

			int hcount = atom.getImplicitHydrogenCount() == CDKConstants.UNSET ? 0 : atom.getImplicitHydrogenCount();
			int charge = atom.getFormalCharge() == CDKConstants.UNSET ? 0 : atom.getFormalCharge();

			double valence;
			// calculate the valence
			valence = ElementCalculation.calculateValence(atom.getSymbol());

			if (bondOrderSum - charge + hcount != valence) {
				double diff = bondOrderSum - charge + hcount;
				if (diff < ElementCalculation.calculateMinimumValence(atom.getSymbol()) ||
						diff > ElementCalculation.calculateMaximumValence(atom.getSymbol())){
					/*System.out.println(ac.getConnectedBondsCount(atom)); 
					 for (IBond b : ac.getConnectedBondsList(atom)) {
					 System.out.println("  "+b.getBegin().getSymbol() + " " +
					 b.getEnd().getSymbol()); } System.out.println(atom.getSymbol() +
					 " bondOrderSum " + bondOrderSum + " charge " + charge + " hcount " + hcount +
					 " valence " + valence + " sum " + ac.getBondOrderSum(atom));*/
					return false;
				}
			}
		}
		return true;
	}

	/**
	 * @param atom
	 * @return
	 */
	public static boolean isValencyOfMoleculeValid(IAtomContainer ac, IAtom atom) {
		if (atom.getSymbol().equals("R")) {
			atom.setImplicitHydrogenCount(0);
			return true;
		}
		if (ac.getConnectedBondsCount(atom) == 0) {
			return true;
		}
		double bondOrderSum = ac.getBondOrderSum(atom);
		int hcount = atom.getImplicitHydrogenCount() == CDKConstants.UNSET ? 0 : atom.getImplicitHydrogenCount();
		int charge = atom.getFormalCharge() == CDKConstants.UNSET ? 0 : atom.getFormalCharge();

		double valence;
		// calculate the valence
		valence = ElementCalculation.calculateValence(atom.getSymbol());

		if (bondOrderSum - charge + hcount != valence) {
			return false;
		}

		return true;
	}

	public static int caclulateHydrogen(IAtomContainer ac, IAtom atom, boolean kekulized) {
		if (atom.getSymbol().equals("R")) {
			atom.setImplicitHydrogenCount(0);
			return -1;
		}
		// double bondOrderSum = ac.getBondOrderSum(atom);
		double bondOrderSum = 0;
		boolean isAromatic = false;
		 /*System.out.println(atom.getSymbol() + " "+ac.getBondOrderSum(atom) + " "
		 +ac.getConnectedBondsCount(atom));*/
		for (IBond bond : ac.getConnectedBondsList(atom)) {
			IBond.Order order = bond.getOrder();
			 //System.out.println(order + " " + order.numeric() + " " +bond.isAromatic());
			if (order.equals(IBond.Order.UNSET) && bond.isAromatic()) {
				bondOrderSum += 1;
				isAromatic = true;
			} else if (!order.equals(IBond.Order.UNSET) && bond.isAromatic()) {
				if (kekulized) {
					bondOrderSum += order.numeric();
				} else {
					bondOrderSum += 1;
					isAromatic = true;
				}
			} else if (!order.equals(IBond.Order.UNSET) && !bond.isAromatic())
				bondOrderSum += order.numeric();
		}
		// adjust aromatic bond score temp fix if c(c)c score should be 4 if cc 3 and 3
		// for cn(C)n (will be adjusted below)
		if (isAromatic)
			bondOrderSum += 1;
		int hcount = atom.getImplicitHydrogenCount() == CDKConstants.UNSET ? 0 : atom.getImplicitHydrogenCount();
		int charge = atom.getFormalCharge() == CDKConstants.UNSET ? 0 : atom.getFormalCharge();

		double valence;
		// calculate the valence
		valence = ElementCalculation.calculateValence(atom.getSymbol());

		if (bondOrderSum - charge + hcount != valence) {
			int diff = (int) (valence - (bondOrderSum -charge + hcount));
			if (diff < 0 && isAromatic) diff += 1;
			int implicitHydrogen2 = atom.getImplicitHydrogenCount() + diff;
			int implicitHydrogen = implicitHydrogen2 > 0 ? implicitHydrogen2 : 0;
			 /*System.out.println("CACLULATE " + atom.getSymbol() + " diff " + diff + 
			 "bondOrderSum " + bondOrderSum + " charge " + charge + " hcount " + hcount + 
			 "valence " + valence + " returned h " + implicitHydrogen);*/
			 
			return implicitHydrogen;
		}
		return -1;
	}

	/**
	 * @param reaction
	 * @param aam
	 * @return
	 * @throws CDKException
	 */
	public static String makeSmiles(IReaction reaction, boolean aam) throws CDKException {
		SmilesGenerator sg;
		if (aam == false)
			sg = new SmilesGenerator(SmiFlavor.Stereo | SmiFlavor.UseAromaticSymbols);
		else
			sg = new SmilesGenerator(SmiFlavor.Stereo | SmiFlavor.AtomAtomMap | SmiFlavor.UseAromaticSymbols);
		//return sg.create(reaction);
		String smi = "";
		for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
			smi+=tools.makeSmiles(ac, false, false, false) + ".";

		}
		smi = smi.substring(0, smi.length() - 1);
		smi+=">";
		for (IAtomContainer ac : reaction.getAgents().atomContainers()) {
			smi+=tools.makeSmiles(ac, false, false, false) + ".";

		}
		smi = smi.substring(0, smi.length() - 1);
		smi+=">";
		for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
			smi+=tools.makeSmiles(ac, false, false, false) + ".";

		}
		smi = smi.substring(0, smi.length() - 1);
		return smi;
	}

	/**
	 * @param ac
	 * @param stereo
	 * @param aam
	 * @param missingHydrogen
	 * @return
	 * @throws CDKException
	 */
	public static String makeSmiles(IAtomContainer ac, boolean stereo, boolean aam, boolean missingHydrogen)
			throws CDKException {
		SmilesGenerator sg;

		if (missingHydrogen)
			addMissingHydrogen(ac, -1, false);

		if (stereo == false && aam == false)
			sg = new SmilesGenerator(SmiFlavor.Canonical);
		else if (stereo == true && aam == false)
			sg = new SmilesGenerator(SmiFlavor.Absolute);
		else if (stereo == false && aam == true)
			sg = new SmilesGenerator(SmiFlavor.AtomAtomMap);
		else
			sg = new SmilesGenerator(SmiFlavor.Absolute | SmiFlavor.AtomAtomMap);

		return sg.create(ac);
	}
}
