package com.nih.tools;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import com.nih.reaction.additionalConstants;

public class HydrogenAdder {

	public static void addMissingHydrogen(IAtomContainer ac, IAtom atom, int maxDepth, boolean kekulized) {
		if (maxDepth > -1) {
			if ((int) atom.getProperty("depth") == maxDepth)
				return;
		}
		if (atom.getSymbol().equals("R")) {
			atom.setImplicitHydrogenCount(0);
			return;
		}

		//double bondOrderSum = ac.getBondOrderSum(atom);
		int bondOrderSum = 0;
		boolean aromaticAdjustment = false;
		for (IBond bond : ac.getConnectedBondsList(atom)) {
			IBond.Order order = bond.getOrder();
			if (order.equals(IBond.Order.UNSET) && bond.isAromatic()){
				bondOrderSum += 1;
				aromaticAdjustment = true;
			}
			else if (!order.equals(IBond.Order.UNSET) && bond.isAromatic()){
				if (kekulized) {
					bondOrderSum += order.numeric();
				}
				else {
					bondOrderSum += 1;
					aromaticAdjustment = true;
				}
			}
			else if (!order.equals(IBond.Order.UNSET) && !bond.isAromatic())
				bondOrderSum += order.numeric();
		}
		//adjust aromatic bond score temp fix if c(c)c score should be 4 if cc 3 and 3 for cn(C)n (will be adjusted below)
		if (aromaticAdjustment) bondOrderSum += 1;

		int charge = atom.getFormalCharge() == CDKConstants.UNSET ? 0 : atom.getFormalCharge();
		int radical = atom.getProperty(additionalConstants.RADICAL) == null ? 0 : atom.getProperty(additionalConstants.RADICAL);

		int valence = (int) (-charge + radical + bondOrderSum);

		Element elem = Element.valueOfIgnoreCase(atom.getSymbol()) ;
		if (atom.isAromatic())
			atom.setImplicitHydrogenCount(implicitHydrogenCount(elem, valence));
		else
			atom.setImplicitHydrogenCount(implicitAromHydrogenCount(elem, valence));
	}
	
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
			int bondOrderSum = 0;
			boolean aromaticAdjustment = false;
			for (IBond bond : ac.getConnectedBondsList(atom)) {
				IBond.Order order = bond.getOrder();
				if (order.equals(IBond.Order.UNSET) && bond.isAromatic()){
					bondOrderSum += 1;
					aromaticAdjustment = true;
				}
				else if (!order.equals(IBond.Order.UNSET) && bond.isAromatic()){
					if (kekulized) {
						bondOrderSum += order.numeric();
					}
					else {
						bondOrderSum += 1;
						aromaticAdjustment = true;
					}
				}
				else if (!order.equals(IBond.Order.UNSET) && !bond.isAromatic())
					bondOrderSum += order.numeric();
			}
			//adjust aromatic bond score temp fix if c(c)c score should be 4 if cc 3 and 3 for cn(C)n (will be adjusted below)
			if (aromaticAdjustment) bondOrderSum += 1;
			
			int charge = atom.getFormalCharge() == CDKConstants.UNSET ? 0 : atom.getFormalCharge();
			int radical = atom.getProperty(additionalConstants.RADICAL) == null ? 0 : atom.getProperty(additionalConstants.RADICAL);

			int valence = (int) (-charge + radical + bondOrderSum);
						
			Element elem = Element.valueOfIgnoreCase(atom.getSymbol()) ;
			if (atom.isAromatic())
				atom.setImplicitHydrogenCount(implicitAromHydrogenCount(elem, valence));
			else
				atom.setImplicitHydrogenCount(implicitHydrogenCount(elem, valence));
		}
	}
	
	/**
     * Determine the implicit hydrogen count of an organic subset atom
     * given its bonded valence. The number of implied hydrogens an 
     * organic (or aromatic) subset atom has is based on it's bonded
     * valence. The valences for the organic elements (B, C, N, O, P,
     * S, F, Cl, Br and I) are defined in the OpenSMILES specification.
     *
     * @param elem Element
     * @param v    bonded valence
     * @return hydrogen count >= 0
     */
    static int implicitHydrogenCount(final Element elem, final int v) {
        switch (elem) {
            case B:
                if (v < 3)  return 3-v;
                break;
            case C:
                if (v < 4)  return 4-v;
                break;
            case N:
            	if (v < 3)  return 3-v;
                break;
            case P:
                if (v <= 3) return 3-v;
                if (v <  5) return 5-v;
                break;
            case O:
                if (v < 2)  return 2-v;
                break;
            case S:
                if (v <= 2) return 2-v;
                if (v <= 4) return 4-v;
                if (v <  6) return 6-v;
                break;
            case Cl:
            	if (v < 1)  return 1;
                break;
            case Br:
            	if (v < 1)  return 1;
                break;
            case I:
            	if (v < 1)  return 1;
                break;
            case F:
                if (v < 1)  return 1;
                break;
		default:
			int val = ElementCalculation.calculateValence(elem)-v;
			return val > -1 ? val : 0;
        }
        return 0;
    }

    /**
     * Determine the implicit hydrogen count of an organic subset atom
     * given its bonded valence. The number of implied hydrogens an 
     * organic (or aromatic) subset atom has is based on it's bonded
     * valence. The valences for the organic elements (B, C, N, O, P,
     * S, F, Cl, Br and I) are defined in the OpenSMILES specification.
     * For aromatic atoms we only check the first level.
     *
     * @param elem Element
     * @param v    bonded valence
     * @return hydrogen count >= 0
     */
    public static int implicitAromHydrogenCount(final Element elem, final int v) {
        switch (elem) {
            case B: // arom?
                if (v < 3) return 3-v;
                break;
            case C:
                if (v < 4) return 4-v;
                break;
            case N:
            	if (v < 3) return 3-v;
                break;
            case P:
                if (v < 3) return 3-v;
                break;
            case O:
                if (v < 2) return 2-v;
                break;
            case S:
                if (v < 2) return 2-v;
                break;
        }
        return 0;
    }
    
  
	
}
