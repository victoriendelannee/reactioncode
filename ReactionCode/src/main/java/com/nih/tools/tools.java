package com.nih.tools;

import static com.nih.reaction.additionalConstants.BOND_CHANGE_INFORMATION;
import static com.nih.reaction.additionalConstants.BOND_CLEAVED;
import static com.nih.reaction.additionalConstants.BOND_MADE;
import static com.nih.reaction.additionalConstants.BOND_ORDER;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.SingleElectron;
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
import org.openscience.cdk.interfaces.IElectronContainer;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.ISingleElectron;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality;
//import org.openscience.cdk.interfaces.IRingSet;
//import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.stereo.Atropisomeric;
import org.openscience.cdk.stereo.DoubleBondStereochemistry;
import org.openscience.cdk.stereo.ExtendedCisTrans;
import org.openscience.cdk.stereo.ExtendedTetrahedral;
import org.openscience.cdk.stereo.Octahedral;
import org.openscience.cdk.stereo.SquarePlanar;
import org.openscience.cdk.stereo.TetrahedralChirality;
import org.openscience.cdk.stereo.TrigonalBipyramidal;

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
		if (!isAromaticityAlreadyPerceived(ac)) {
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
	}
		

		
	
	private static boolean isAromaticityAlreadyPerceived(IAtomContainer ac) {
		for (IBond bond : ac.bonds()) {
			if (bond.isAromatic())
				return true;
		}
		return false;
	}


	public static boolean checkProductValidity(IAtomContainerSet acset) {
		// check if generated product is valid (only check modified atoms)
		for (IAtomContainer ac : acset.atomContainers()) {
			boolean isValencyValid = isValencyOfMoleculeValid(ac, true);
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
	public static boolean isValencyOfMoleculeValid(IAtomContainer ac, boolean kekulize) {
		boolean mustReturn = false;
		if (kekulize) {
			try {
				Kekulization.kekulize(ac);
				if (!Kekulization2.kekulize(ac))
					return false;
			} catch (CDKException e) {
				// TODO Auto-generated catch block
				// e.printStackTrace();
				mustReturn = true;
				// return true;
			}
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
			int radical = atom.getProperty(additionalConstants.RADICAL) == null ? 0 : atom.getProperty(additionalConstants.RADICAL);


			// calculate the valence
			double cVal = bondOrderSum - charge + hcount + radical;
			return checkValenceValidity(Element.valueOfIgnoreCase(atom.getSymbol()) , cVal) ;
			//OLD Method
			/*
			double valence;
			//calculate the valence
			valence = ElementCalculation.calculateValence(atom.getSymbol());

			if (bondOrderSum - charge + hcount + radical != valence) {
				double diff = bondOrderSum - charge + hcount + radical;
				if (diff < ElementCalculation.calculateMinimumValence(atom.getSymbol()) ||
						diff > ElementCalculation.calculateMaximumValence(atom.getSymbol())){
					return false;
				}
			}
			*/
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
		int radical = atom.getProperty(additionalConstants.RADICAL) == null ? 0 : atom.getProperty(additionalConstants.RADICAL);

		
		// calculate the valence
		double cVal = bondOrderSum - charge + hcount + radical;
		return checkValenceValidity(Element.valueOfIgnoreCase(atom.getSymbol()) , cVal) ;
		//OLD METHOD
		/*
		 double valence = ElementCalculation.calculateValence(atom.getSymbol());
		if (cVal != valence) {
			
			return false;
		}

		return true;
		*/
	}

	public static void caclulateHydrogen(IAtomContainer ac, IAtom atom, boolean kekulized) {
		HydrogenAdder.addMissingHydrogen(ac, atom, -1, kekulized);
	}
	
	
	
	/**
	 * Deduce radicals
	 * @param ac
	 * @param kekulized
	 */
	public static void radicalize(IAtomContainer ac, boolean kekulized) {
		for (IAtom atom : ac.atoms()) {
			int radical = radicalize(ac, atom, kekulized);
			atom.setProperty(additionalConstants.RADICAL, radical);
		}
	}
	
	/**
	 * Deduce radicals
	 * @param ac
	 * @param atom
	 * @param kekulized
	 * @return
	 */
	public static int radicalize(IAtomContainer ac, IAtom atom, boolean kekulized) {
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
			int diff = (int) (valence - (bondOrderSum - charge + hcount));
			if (diff < 0 && isAromatic) diff += 1;
			int radicalCount2 = atom.getImplicitHydrogenCount() + diff;
			int radicalCount = radicalCount2 > 0 ? radicalCount2 : 0;
			 /*System.out.println("CACLULATE " + atom.getSymbol() + " diff " + diff + 
			 "bondOrderSum " + bondOrderSum + " charge " + charge + " hcount " + hcount + 
			 "valence " + valence + " returned radical " + radicalCount);*/
			if (radicalCount > 0) {
                for (int j = 0; j < radicalCount; j++) {
                    ac.addSingleElectron(ac.getBuilder().newInstance(ISingleElectron.class,
                    		atom));
                }
                return radicalCount;
			}
			 
			return -1;
		}
		return -1;
	}
	

    private static boolean checkValenceValidity(final Element elem, final double cVal) {
        switch (elem) {
            case P:
                if (cVal <= 3) return 3 == cVal;
                if (cVal >  3) return 5 == cVal;
                break;
            case S:
                if (cVal <= 2) return 2 == cVal;
                if (cVal <= 4) return 4 == cVal;
                if (cVal > 4) return 6 == cVal;
                break;
		default:
			return ElementCalculation.calculateValence(elem) == cVal;
        }
		return false;
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
			smi+=tools.makeSmiles(ac, true, aam, false) + ".";

		}
		smi = smi.substring(0, smi.length() - 1);
		smi+=">";
		for (IAtomContainer ac : reaction.getAgents().atomContainers()) {
			smi+=tools.makeSmiles(ac, true, aam, false) + ".";

		}
		if (reaction.getAgents().getAtomContainerCount() > 0)
			smi = smi.substring(0, smi.length() - 1);
		smi+=">";
		for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
			smi+=tools.makeSmiles(ac, true, aam, false) + ".";

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
			HydrogenAdder.addMissingHydrogen(ac, -1, false);

		if (stereo == false && aam == false)
			sg = new SmilesGenerator(SmiFlavor.Canonical | SmiFlavor.UseAromaticSymbols);
		else if (stereo == true && aam == false)
			sg = new SmilesGenerator(SmiFlavor.Absolute | SmiFlavor.UseAromaticSymbols);
		else if (stereo == false && aam == true)
			sg = new SmilesGenerator(SmiFlavor.AtomAtomMap | SmiFlavor.UseAromaticSymbols);
		else
			sg = new SmilesGenerator(SmiFlavor.Absolute | SmiFlavor.AtomAtomMap | SmiFlavor.UseAromaticSymbols);

		return sg.create(ac);
	}
	
	public static IReaction convertReactionContainingIQueryAtomContainer(IReaction reaction) {
		Set<IAtom> reactionCenter = new HashSet<IAtom>();
		Set<IBond> bondsCleaved = new HashSet<IBond>();
		Set<IBond> bondsFormed = new HashSet<IBond>();
		Set<IBond> bondsOrder = new HashSet<IBond>();
		
		IReaction newReaction = DefaultChemObjectBuilder.getInstance().newInstance(IReaction.class);
		IAtomContainerSet newReactants = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		IAtomContainerSet newProducts = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		//old new
		for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
			IAtomContainer nac = convertIQueryAtomContainer(ac);
			reactionCenter.addAll(nac.getProperty("reactionCenter"));
			bondsCleaved.addAll(nac.getProperty("bondsCleaved"));
			bondsFormed.addAll(nac.getProperty("bondsFormed"));
			bondsOrder.addAll(nac.getProperty("bondsOrder"));
			newReactants.addAtomContainer(nac);
		}
		for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
			IAtomContainer nac = convertIQueryAtomContainer(ac);
			reactionCenter.addAll(nac.getProperty("reactionCenter"));
			bondsCleaved.addAll(nac.getProperty("bondsCleaved"));
			bondsFormed.addAll(nac.getProperty("bondsFormed"));
			bondsOrder.addAll(nac.getProperty("bondsOrder"));
			newProducts.addAtomContainer(nac);
		}
		newReaction.setReactants(newReactants);
		newReaction.setProducts(newProducts);
		
		newReaction.setID(reaction.getID());
		newReaction.setProperty("reactionCenter", reactionCenter);
		newReaction.setProperty("bondsCleaved", bondsCleaved);
		newReaction.setProperty("bondsFormed", bondsFormed);
		newReaction.setProperty("bondsOrder", bondsOrder);
		return newReaction;
	}
	
	private static IAtomContainer convertIQueryAtomContainer(IAtomContainer ac) {
		Set<IAtom> reactionCenter = new HashSet<IAtom>();
		Set<IBond> bondsCleaved = new HashSet<IBond>();
		Set<IBond> bondsFormed = new HashSet<IBond>();
		Set<IBond> bondsOrder = new HashSet<IBond>();
		
		//old new
		Map<IAtom,IAtom> index1;
		Map<IBond,IBond> index2;
		index1 = new HashMap<IAtom,IAtom>();
		index2 = new HashMap<IBond,IBond>();
		IAtomContainer nac = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
		for (IAtom a : ac.atoms()) {
			IAtom na = new Atom(a.getSymbol());
			na.setFormalCharge(a.getFormalCharge());
			na.setIsAromatic(a.isAromatic());
			na.setIsInRing(a.isInRing());
			na.setID(a.getID());
			na.setProperties(a.getProperties());
			na.setImplicitHydrogenCount(a.getImplicitHydrogenCount());
			nac.addAtom(na);
			int singleElecronCount = a.getProperty(additionalConstants.RADICAL) == null ? 0 : a.getProperty(additionalConstants.RADICAL);
			for (int j = 0; j < singleElecronCount; j++) {
				ISingleElectron se = new SingleElectron(na);
				nac.addSingleElectron(se);
			}
			index1.put(a, na);
		}
		for (IBond b : ac.bonds()) {
			IBond nb = new Bond();
			nb.setAtom(index1.get(b.getBegin()), 0);
			nb.setAtom(index1.get(b.getEnd()), 1);
			nb.setIsAromatic(b.isAromatic());
			nb.setIsInRing(b.isInRing());
			nb.setOrder(b.getOrder());
			nb.setStereo(b.getStereo());
			nb.setProperties(b.getProperties());
			nac.addBond(nb);
			index2.put(b, nb);
			if (nb.getProperty(BOND_CHANGE_INFORMATION) != null) {
				if ((int)nb.getProperty(BOND_CHANGE_INFORMATION) == BOND_CLEAVED) {
					bondsCleaved.add(nb);
					reactionCenter.add(nb.getBegin());
					reactionCenter.add(nb.getEnd());
				}
				if ((int)nb.getProperty(BOND_CHANGE_INFORMATION) == BOND_MADE) {
					bondsFormed.add(nb);
					reactionCenter.add(nb.getBegin());
					reactionCenter.add(nb.getEnd());
				}
				if ((int)nb.getProperty(BOND_CHANGE_INFORMATION) == BOND_ORDER) {
					bondsOrder.add(nb);
					reactionCenter.add(nb.getBegin());
					reactionCenter.add(nb.getEnd());
				}
			}
		}
		for (IStereoElement se : ac.stereoElements()) {
			if (se instanceof TetrahedralChirality) {
				IAtom focus = index1.get(se.getFocus());
				IAtom[] old = ((TetrahedralChirality) se).getLigands();
				IAtom[] ligands = new IAtom[old.length];
				for (int i = 0; i < old.length; i++) {
					ligands[i] = index1.get(old[i]);
				}
				nac.addStereoElement(
						new TetrahedralChirality((IAtom)focus, ligands, se.getConfig()));
			}
			else if (se instanceof ExtendedTetrahedral) {
				IAtom focus = index1.get(se.getFocus());
				IAtom[] old = ((TetrahedralChirality) se).getLigands();
				IAtom[] ligands = new IAtom[old.length];
				for (int i = 0; i < old.length; i++) {
					ligands[i] = index1.get(old[i]);
				}
				nac.addStereoElement(
						new ExtendedTetrahedral((IAtom)focus, ligands, se.getConfig()));
			}
			else if (se instanceof SquarePlanar) {
				IAtom focus = index1.get(se.getFocus());
				List<IAtom> old = ((SquarePlanar) se).getCarriers();
				IAtom[] ligands = new IAtom[old.size()];
				for (int i = 0; i < old.size(); i++) {
					ligands[i] = index1.get(old.get(i));
				}
				nac.addStereoElement(
						new SquarePlanar((IAtom)focus, ligands, se.getConfig()));
			}
			else if (se instanceof TrigonalBipyramidal) {
				IAtom focus = index1.get(se.getFocus());
				List<IAtom> old = ((TrigonalBipyramidal) se).getCarriers();
				IAtom[] ligands = new IAtom[old.size()];
				for (int i = 0; i < old.size(); i++) {
					ligands[i] = index1.get(old.get(i));
				}
				nac.addStereoElement(
						new TrigonalBipyramidal((IAtom)focus, ligands, se.getConfig()));
			}
			else if (se instanceof Octahedral) {
				IAtom focus = index1.get(se.getFocus());
				List<IAtom> old = ((Octahedral) se).getCarriers();
				IAtom[] ligands = new IAtom[old.size()];
				for (int i = 0; i < old.size(); i++) {
					ligands[i] = index1.get(old.get(i));
				}
				nac.addStereoElement(
						new Octahedral((IAtom)focus, ligands, se.getConfig()));
			}
			else if (se instanceof DoubleBondStereochemistry) {
				IBond focus = index2.get(se.getFocus());
				IBond[] old = ((DoubleBondStereochemistry) se).getBonds();
				IBond[] ligands = new IBond[old.length];
				for (int i = 0; i < old.length; i++) {
					ligands[i] = index2.get(old[i]);
				}
				nac.addStereoElement(
						new DoubleBondStereochemistry((IBond)focus, ligands, se.getConfig()));
			}
			else if (se instanceof ExtendedCisTrans) {
				IBond focus = index2.get(se.getFocus());
				List<IBond> old = ((ExtendedCisTrans) se).getCarriers();
				IBond[] ligands = new IBond[old.size()];
				for (int i = 0; i < old.size(); i++) {
					ligands[i] = index2.get(old.get(i));
				}
				nac.addStereoElement(
						new ExtendedCisTrans((IBond)focus, ligands, se.getConfig()));
			}
			else if (se instanceof Atropisomeric) {
				IBond focus = index2.get(se.getFocus());
				List<IAtom> old = ((Atropisomeric) se).getCarriers();
				IAtom[] ligands = new IAtom[old.size()];
				for (int i = 0; i < old.size(); i++) {
					ligands[i] = index1.get(old.get(i));
				}
				nac.addStereoElement(
						new Atropisomeric((IBond)focus, ligands, se.getConfig()));
			}
		}
		nac.setProperty("reactionCenter", reactionCenter);
		nac.setProperty("bondsCleaved", bondsCleaved);
		nac.setProperty("bondsFormed", bondsFormed);
		nac.setProperty("bondsOrder", bondsOrder);
		return nac;
	}
	
	//old calculation method
	public static int caclulateHydrogen2(IAtomContainer ac, IAtom atom, boolean kekulized) {
		if (atom.getSymbol().equals("R")) {
			atom.setImplicitHydrogenCount(0);
			return -1;
		}

		double bondOrderSum = 0;
		boolean isAromatic = false;

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
		int radical = atom.getProperty(additionalConstants.RADICAL) == null ? 0 : atom.getProperty(additionalConstants.RADICAL);

		double valence;
		// calculate the valence
		valence = ElementCalculation.calculateValence(atom.getSymbol());

		if (bondOrderSum - charge + hcount + radical != valence) {
			int diff = (int) (valence - (bondOrderSum - charge + hcount + radical));
			if (diff < 0 && isAromatic) diff += 1;
			int implicitHydrogen2 = atom.getImplicitHydrogenCount() + diff;
			int implicitHydrogen = implicitHydrogen2 > 0 ? implicitHydrogen2 : 0;

			return implicitHydrogen;
		}
		return -1;
	}
	
	//old calculation method
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
			double bondOrderSum = 0;
			boolean isAromatic = false;

			for (IBond bond : ac.getConnectedBondsList(atom)) {
				IBond.Order order = bond.getOrder();
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
			int radical = atom.getProperty(additionalConstants.RADICAL) == null ? 0 : atom.getProperty(additionalConstants.RADICAL);

			double valence;
			// calculate the valence
			valence = ElementCalculation.calculateValence(atom.getSymbol());

			if (bondOrderSum - charge + hcount +radical != valence) {
				int diff = (int) (valence - (bondOrderSum - charge + hcount + radical));
				if (diff < 0 && isAromatic) diff += 1;
				int implicitHydrogen2 = atom.getImplicitHydrogenCount() + diff;
				int implicitHydrogen = implicitHydrogen2 > 0 ? implicitHydrogen2 : 0;
				atom.setImplicitHydrogenCount(implicitHydrogen);
			}
			else {
				if (atom.getImplicitHydrogenCount() == null)
					atom.setImplicitHydrogenCount(0);
			}
		}
	}
}
