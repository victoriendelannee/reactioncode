package com.nih.tools;

import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.CDKConstants;

/*
 * Copyright (c) 2014 European Bioinformatics Institute (EMBL-EBI)
 *                    John May <jwmay@users.sf.net>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version. All we ask is that proper credit is given
 * for our work, which includes - but is not limited to - adding the above
 * copyright notice to the beginning of your source code files, and to any
 * copyright notice that you may distribute with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 U
 */

import static com.nih.reaction.additionalConstants.BOND_CHANGE_INFORMATION;

import org.openscience.cdk.config.Elements;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.Intractable;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.Matching2;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.RingPartitioner;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import com.nih.reaction.additionalConstants;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import static org.openscience.cdk.CDKConstants.ISAROMATIC;
import static org.openscience.cdk.graph.GraphUtil.EdgeToBondMap;
import static org.openscience.cdk.interfaces.IBond.Order.DOUBLE;
import static org.openscience.cdk.interfaces.IBond.Order.SINGLE;
import static org.openscience.cdk.interfaces.IBond.Order.UNSET;

/**
 * Assign a Kekulé representation to the aromatic systems of a compound. Input
 * from some file-formats provides some bonds as aromatic / delocalised bond
 * types. This method localises the electrons and assigns single and double
 * bonds. Different atom and bond orderings may produce distinct but valid
 * Kekulé forms. Only bond orders are adjusted and any aromatic flags will
 * remain untouched.
 * 
 *
 * The procedure requires that all atoms have defined implicit hydrogens counts
 * and formal charges. If this information is not present it should be assigned
 * first. 
 *
 * For some inputs it may not be possible to assign a Kekulé form. In general
 * theses cases are rare but usually occur for one of two reasons.
 * 1) Missing / ambiguous implicit hydrogens, this is fundamental to determining the
 * Kekulé form and if guessed may be wrong. Some formats (e.g. molfile) can not
 * include the exact number of implicit hydrogens attached to atom whilst others
 * may omit it or optionally skip encoding. The typical example is found in the
 * example for 1H-pyrrole, a correct SMILES encoding should include the hydrogen
 * on the aromatic nitrogen '[nH]1cccc1' (not: 'n1cccc1').
 * 2) The aromaticity perception algorithm has allowed atoms with abnormal
 * valence. This usually happens when a non-convalent bond has be <i>upgraded</i>
 * to a sigma bond during format conversion. 
 *
 * @author John May
 * @cdk.keyword kekule
 * @cdk.keyword kekulize
 * @cdk.keyword dearomatize
 * @cdk.keyword aromatic
 * @cdk.keyword fix bond orders
 * @cdk.keyword deduce bond orders
 */
public final class Kekulization2 {

	public static boolean kekulizeUsingReference(final IAtomContainer reference, final IAtomContainer target) throws CDKException {
		Map<Integer,IAtom> index = new HashMap<Integer,IAtom>();
		for (IAtom atom : target.atoms()) {
			index.put(atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING), atom);
		}
		
		for (IBond rBond : reference.bonds()) {
			if (rBond.isAromatic() && rBond.getProperty(BOND_CHANGE_INFORMATION) == null) {
				int aamB = rBond.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING);
				int aamE = rBond.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING);
				IBond tBond = target.getBond(index.get(aamB), index.get(aamE));
				if (tBond != null) {
					if (!tBond.isAromatic())
						continue;
					if ((boolean)tBond.getBegin().getProperty(additionalConstants.REACTION_CENTER) == false &&
							(boolean)tBond.getEnd().getProperty(additionalConstants.REACTION_CENTER) == false) {
						tBond.getBegin().setImplicitHydrogenCount(rBond.getBegin().getImplicitHydrogenCount());
						tBond.getEnd().setImplicitHydrogenCount(rBond.getEnd().getImplicitHydrogenCount());
						tBond.setOrder(rBond.getOrder());
					}
				}
			}
		}
		return kekulize(target);
	}
    /**
     * Assign a Kekulé representation to the aromatic systems of a compound.
     *
     * @param ac structural representation
     * @throws CDKException a Kekulé form could not be assigned
     */
    public static boolean kekulize(final IAtomContainer ac) throws CDKException {
    	//check if all aromatic are rings are closed
    	if (!areAllAromaticRingsClosed(ac))
    		return false;
        // storage of pairs of atoms that have pi-bonded
        final Matching2 matching = Matching2.withCapacity(ac.getAtomCount());

        // extract data structures for efficient access
        final IAtom[] atoms = AtomContainerManipulator.getAtomArray(ac);
        final EdgeToBondMap bonds = EdgeToBondMap.withSpaceFor(ac);
        final int[][] graph = GraphUtil.toAdjList(ac, bonds);

        // determine which atoms are available to have a pi bond placed
        final BitSet available = available(graph, atoms, bonds);

        // attempt to find a perfect matching such that a pi bond is placed
        // next to each available atom. if not found the solution is ambiguous
        if (!matching.perfect(graph, available)) {
        	return kekulize2(ac);
        	//throw new CDKException("Cannot assign Kekulé structure without randomly creating radicals.");
        }
        // propegate bond order information from the matching
        for (final IBond bond : ac.bonds()) {
            if (bond.getOrder() == UNSET && bond.isAromatic()) bond.setOrder(SINGLE);
        }
        for (int v = available.nextSetBit(0); v >= 0; v = available.nextSetBit(v + 1)) {
            final int w = matching.other(v);
            final IBond bond = bonds.get(v, w);

            // sanity check, something wrong if this happens
            if (bond.getOrder().numeric() > 1)
                throw new CDKException(
                        "Cannot assign Kekulé structure, non-sigma bond order has already been assigned?");

            bond.setOrder(IBond.Order.DOUBLE);
            available.clear(w);
        }
        return true;
    }
    
    /**
     * Try other methods if original CDK kekulization methof failed
     * @param ac
     * @return
     * @throws CDKException
     */
    public static boolean kekulize2(final IAtomContainer ac) throws CDKException {
    	//randomly set one heteroatom
    	if (!saturateOneRandomlySelectedHeteroAtom(ac))
    		return kekulizeByRingSystem(ac);
    	 // storage of pairs of atoms that have pi-bonded
        final Matching2 matching = Matching2.withCapacity(ac.getAtomCount());

        // extract data structures for efficient access
        final IAtom[] atoms = AtomContainerManipulator.getAtomArray(ac);
        final EdgeToBondMap bonds = EdgeToBondMap.withSpaceFor(ac);
        final int[][] graph = GraphUtil.toAdjList(ac, bonds);

        // determine which atoms are available to have a pi bond placed
        final BitSet available = available(graph, atoms, bonds);

        // attempt to find a perfect matching such that a pi bond is placed
        // next to each available atom. if not found the solution is ambiguous
        if (!matching.perfect(graph, available)) {
        	return false;
        	//throw new CDKException("Cannot assign Kekulé structure without randomly creating radicals.");
        }
        // propegate bond order information from the matching
        for (final IBond bond : ac.bonds()) {
            if (bond.getOrder() == UNSET && bond.isAromatic()) bond.setOrder(SINGLE);
        }
        for (int v = available.nextSetBit(0); v >= 0; v = available.nextSetBit(v + 1)) {
            final int w = matching.other(v);
            final IBond bond = bonds.get(v, w);

            // sanity check, something wrong if this happens
            if (bond.getOrder().numeric() > 1)
                throw new CDKException(
                        "Cannot assign Kekulé structure, non-sigma bond order has already been assigned?");

            bond.setOrder(IBond.Order.DOUBLE);
            available.clear(w);
        }
        return true;
        
    	//return forceKekulize(ac);
    }
    
    /**
	 * Split the molecule by ring system and kekulized according to each ring system
	 * @param ac
	 * @return
	 * @throws CDKException
	 */
	public static boolean kekulizeByRingSystem(final IAtomContainer ac) throws CDKException {
		CycleFinder cf = Cycles.mcb();
    	Cycles   cycles = cf.find(ac);
		IRingSet rings  = cycles.toRingSet();
		List<IRingSet> systems = RingPartitioner.partitionRings(rings);
		boolean fullyKekulzed = true;
		for (IRingSet systemSet : systems) {
			Set<IAtom> toRemove1 = new HashSet<IAtom>();
			Set<IBond> toRemove2 = new HashSet<IBond>();
			//MARCHE PAS soit partitionner soit faire propre algo
			IAtomContainer system  = RingPartitioner.convertToAtomContainer(systemSet);
			boolean isAromatic = false;
			for (IAtom a : system.atoms()) {
				//to prevent the same systym being computed multiple times 
				if (a.isAromatic()) {
					isAromatic = true;
					a.setProperty("bondOrderSum2", (int)ac.getBondOrderSum(a));
					for (IBond bond : ac.getConnectedBondsList(a)) {
						if (!system.contains(bond) && !bond.isInRing()) {
							IBond newBond = new Bond();
							newBond.setOrder(bond.getOrder());
							newBond.setIsAromatic(bond.isAromatic());
							newBond.setIsInRing(false);
							IAtom fake = new Atom();
							fake.setIsAromatic(false);
							fake.setIsInRing(false);
							if (!system.contains(bond.getBegin())) {
								fake.setSymbol(bond.getEnd().getSymbol());
								fake.setImplicitHydrogenCount(bond.getEnd().getImplicitHydrogenCount());
								newBond.setAtom(bond.getEnd(), 0);
								newBond.setAtom(fake, 1);
							}
							if (!system.contains(bond.getEnd())) {
								fake.setSymbol(bond.getBegin().getSymbol());
								fake.setImplicitHydrogenCount(bond.getBegin().getImplicitHydrogenCount());
								newBond.setAtom(bond.getBegin(), 0);
								newBond.setAtom(fake, 1);
							}
							toRemove1.add(fake);
							toRemove2.add(newBond);
						}
					}
				}
			}
			if (!isAromatic)
				continue;
			for (IAtom toRemove : toRemove1) {
				system.addAtom(toRemove);
			}
			for (IBond toRemove : toRemove2) {
				system.addBond(toRemove);
			}

			//check if all aromatic are rings are closed
	    	if (!areAllAromaticRingsClosed(system)){
	    		fullyKekulzed = false;
	    		continue;
	    	}
	    	
	    	//case where a system is already valid (made by ref)
			if (checkSolution(system))
				continue;
			
	        // storage of pairs of atoms that have pi-bonded
	        final Matching2 matching = Matching2.withCapacity(system.getAtomCount());

	        // extract data structures for efficient access
	        final IAtom[] atoms = AtomContainerManipulator.getAtomArray(system);
	        final EdgeToBondMap bonds = EdgeToBondMap.withSpaceFor(system);
	        final int[][] graph = GraphUtil.toAdjList(system, bonds);

	        // determine which atoms are available to have a pi bond placed
	        final BitSet available = available(graph, atoms, bonds);

	        // attempt to find a perfect matching such that a pi bond is placed
	        // next to each available atom. if not found the solution is ambiguous
	        if (!matching.perfect(graph, available)) {
	        	if (findSolution(system, available, bonds, matching))
	        		fullyKekulzed = false;
	        	continue;
	        	//throw new CDKException("Cannot assign Kekulé structure without randomly creating radicals.");
	        }

	        // propegate bond order information from the matching
	        for (final IBond bond : system.bonds()) {
	            if (bond.getOrder() == UNSET && bond.isAromatic()) bond.setOrder(SINGLE);
	        }
	        for (int v = available.nextSetBit(0); v >= 0; v = available.nextSetBit(v + 1)) {
	            final int w = matching.other(v);
	            final IBond bond = bonds.get(v, w);
	            // sanity check, something wrong if this happens
	            if (bond.getOrder().numeric() > 1)
	                throw new CDKException(
	                        "Cannot assign Kekulé structure, non-sigma bond order has already been assigned?");

	            bond.setOrder(IBond.Order.DOUBLE);
	            available.clear(w);
	        }
		}
		return fullyKekulzed;
	}
    
    private static boolean findSolution(IAtomContainer ac, BitSet available, EdgeToBondMap bonds, Matching2 matching) {
		//annotate bonds, which can be double bonds
		List<IAtom> availableAtomsL = new ArrayList<IAtom>();
		int best = 99;
		for (int v = available.nextSetBit(0); v >= 0; v = available.nextSetBit(v + 1)) {
			IAtom atom = ac.getAtom(v);
			if (atom.isInRing()){
				int count = ac.getConnectedBondsCount(atom);
				if (count <= best) {
					best = count;
					availableAtomsL.add(0, atom);
				}
				else
					availableAtomsL.add(atom);
			}
        }
		for (IAtom a : ac.atoms()) {
			a.setProperty("kekulizeFlag", false);
		}
		Set<IAtom> availableAtoms = new HashSet<IAtom>(availableAtomsL);
		
		for (IAtom atom : availableAtoms) {
			Set<IAtom> current = new HashSet<IAtom>();
			current.add(atom);
			findSolution(ac, new HashSet<IAtom>(availableAtoms), current, false);
			if (checkSolution(ac))
				return true;
			restoreOriginalBondConf(ac);
			findSolution(ac, new HashSet<IAtom>(availableAtoms), current, true);
			if (checkSolution(ac))
				return true;
			restoreOriginalBondConf(ac);
		}
		return false;
		
	}
	
	/**
	 * restore the original bond order
	 * @param ac
	 */
	private static void restoreOriginalBondConf(IAtomContainer ac) {
		for (IBond bond : ac.bonds()) {
			if (bond.getOrder().equals(IBond.Order.DOUBLE) && bond.isAromatic() && bond.isInRing()){
				bond.setOrder(IBond.Order.SINGLE);
			}
		}
	}
	
	/**
	 * Check if the generated solution is valid
	 * N.planar3: 2 electrons
		N.minus.planar3: 2 electrons
		N.amide: 2 electrons
		S.2: 2 electrons
		S.planar3: 2 electrons
		C.minus.planar: 2 electrons
		O.planar3: 2 electrons
		N.sp2.3: 1 electron
		C.sp2: 1 electron
	 * @param ac
	 * @return
	 */
	private static boolean checkSolution(IAtomContainer ac) {
		int e = 0;
		for (IAtom atom : ac.atoms()) {
			if (!atom.getSymbol().equals("H") && atom.isAromatic()) {
				int charge = atom.getFormalCharge() == CDKConstants.UNSET ? 0 : atom.getFormalCharge();
				int radical = atom.getProperty(additionalConstants.RADICAL) == null ? 0 : atom.getProperty(additionalConstants.RADICAL);
				
				if (atom.getSymbol().equals("C") && atom.getFormalCharge() == 0)
					continue;
				
				int bondOrderSum = 0;
				for (IBond bond : ac.getConnectedBondsList(atom)) {
					IBond.Order order = bond.getOrder();
					bondOrderSum += order.numeric();
				}

				int valence = (int) (-charge + radical + bondOrderSum);
				
				Element elem = Element.valueOfIgnoreCase(atom.getSymbol()) ;

				e += (HydrogenAdder.implicitAromHydrogenCount(elem, valence)*2);
				
				//S 2 electrons and O 2 electrons
				if (atom.getSymbol().equals("O") || atom.getSymbol().equals("S"))
					e+=2;
			}
			atom.setProperty("kekulizeFlag", false);
		}
		for (IBond bond : ac.bonds()) {
			if (bond.getOrder().equals(IBond.Order.DOUBLE) && bond.isAromatic() && bond.isInRing()){
				e+=2;
			}
		}
		return (e - 2) % 4 == 0;
	}
	
	/**
	 * find a valid kekule
	 * @param ac
	 * @param availableAtoms
	 * @param current
	 * @param first
	 * @return
	 */
	private static boolean findSolution(IAtomContainer ac, Set<IAtom> availableAtoms, Set<IAtom> current, boolean first) {
		Set<IAtom> next = new HashSet<IAtom>();
		for (IAtom atom : current) {
			atom.setProperty("kekulizeFlag", true);
			boolean isAttributed = false;
			if (!availableAtoms.contains(atom))
				isAttributed = true;
			availableAtoms.remove(atom);
			//try other direction
			if (first) 
				isAttributed = true;
			List<IBond> bonds = ac.getConnectedBondsList(atom);
			for (int i = 0; i < bonds.size(); i++) {
				IBond bond = bonds.get(i);
				IAtom other = bond.getOther(atom);
				if (!other.isInRing() || (boolean)other.getProperty("kekulizeFlag") == true)
					continue;
				else 
					next.add(other);
				if (first && i == 0)
					isAttributed = true;
				else if (first && i == 1)
					isAttributed = false;
				if (availableAtoms.contains(other)) {
					if (!isAttributed) {
						bond.setOrder(IBond.Order.DOUBLE);
						isAttributed = true;
						availableAtoms.remove(other);
					}
				}
			}
		}
		if (next.isEmpty())
			return true;
		else
			return findSolution(ac, availableAtoms, next, false);
	}
	
    
    /**
     * Check if all rings are closed
     * @param ac
     * @return
     */
    private static boolean areAllAromaticRingsClosed(IAtomContainer ac) {
    	Cycles.markRingAtomsAndBonds(ac);
    	for (IBond bond : ac.bonds()) {
    		if (bond.isAromatic() && !bond.isInRing())
    			return false;
    	}
    	return true;
    }
    
    /**
     * Saturate one hetero atom in order to kekulized
     * @param ac
     * @return
     */
    private static boolean saturateOneRandomlySelectedHeteroAtom(IAtomContainer ac) {
    	CycleFinder cf = Cycles.mcb();
    	try {
			Cycles cycles = cf.find(ac);
			IRingSet rings = cycles.toRingSet();
			Map<String, Set<IAtom>> candidates = new HashMap<String, Set<IAtom>>();

			for (IAtom atom : ac.atoms()) {
				if (atom.isAromatic() && !isSaturate(ac, atom)) {
					String symbol = atom.getSymbol();
					if (!symbol.equals("C") && !symbol.equals("H")) {
						if (candidates.containsKey(symbol)) {
							Set<IAtom> set = candidates.get(symbol);
							set.add(atom);
							candidates.put(symbol, set);
						} else {
							Set<IAtom> set = new HashSet<IAtom>();
							set.add(atom);
							candidates.put(symbol, set);
						}
					}
				}
			}
			Map<String, List<IAtom>> sortedCandidates = sortMapByRing(ac, rings, candidates);
			if (sortedCandidates.containsKey("N")) {
				for (IAtom atom : sortedCandidates.get("N")) {
					HydrogenAdder.addMissingHydrogen(ac, atom, -1, true);
					if (!testMatching(ac))
						atom.setImplicitHydrogenCount(0);
					else
						return true;
				}
			}
			if (sortedCandidates.containsKey("B")) {
				for (IAtom atom : sortedCandidates.get("B")) {
					HydrogenAdder.addMissingHydrogen(ac, atom, -1, true);
					if (!testMatching(ac))
						atom.setImplicitHydrogenCount(0);
					else
						return true;
				}
			}
			if (sortedCandidates.containsKey("P")) {
				for (IAtom atom : sortedCandidates.get("P")) {
					HydrogenAdder.addMissingHydrogen(ac, atom, -1, true);
					if (!testMatching(ac))
						atom.setImplicitHydrogenCount(0);
					else
						return true;
				}
			}
    	}catch (Intractable e) {
    		// ignore error - MCB should never be intractable
    		return false;
    	}
    	return false;
    }
    
    /**
     * sort atom in list by comparing the ring size (atom in smaller rings first) and ReactionCenter info
     * @param ac
     * @param rings
     * @param candidates
     * @return
     */
    private static Map<String,List<IAtom>> sortMapByRing(IAtomContainer ac, IRingSet rings, Map<String,Set<IAtom>> candidates) {
    	Map<String,List<IAtom>> newCandidates = new HashMap<String,List<IAtom>>();
    	for (Entry<String,Set<IAtom>> e : candidates.entrySet()) {
    		String symbol = e.getKey();
    		List<IAtom> atoms = new ArrayList<IAtom>();
    		for (IAtom atom : e.getValue()) {
    			int biggest = 0;
    			for (IAtomContainer ring : rings.getRings(atom).atomContainers()) {
    				int count = ring.getAtomCount();
    				if (count > biggest)
    					biggest = count;
    			}
    			///highest priority
    			boolean isReactionCenter = false;
    			if ((boolean)atom.getProperty(additionalConstants.REACTION_CENTER) == true)
    				isReactionCenter = true;
    			atom.setProperty("biggestRing", biggest);
    			int i;
    			for (i = 0; i < atoms.size(); i++) {
    				IAtom ranked = atoms.get(i);
    				boolean isRc = (boolean)ranked.getProperty(additionalConstants.REACTION_CENTER);
    				if (isReactionCenter) {
    					if (isRc == false) {
    						break;
    					}
    				}
    				if (biggest < (int)ranked.getProperty("biggestRing")) {
    					if (isRc && !isReactionCenter) 
    						continue;
    					break;
    				}
    			}
    			atoms.add(i, atom);
     		}
    		newCandidates.put(symbol, atoms);
    	}
		return newCandidates;
    }
    
    private static boolean testMatching(IAtomContainer ac) {
    	final Matching2 matching = Matching2.withCapacity(ac.getAtomCount());

        // extract data structures for efficient access
        final IAtom[] atoms = AtomContainerManipulator.getAtomArray(ac);
        final EdgeToBondMap bonds = EdgeToBondMap.withSpaceFor(ac);
        final int[][] graph = GraphUtil.toAdjList(ac, bonds);

        // determine which atoms are available to have a pi bond placed
        final BitSet available = available(graph, atoms, bonds);

        // attempt to find a perfect matching such that a pi bond is placed
        // next to each available atom. if not found the solution is ambiguous
        if (!matching.perfect(graph, available)) 
        	return false;
        else 
        	return true;
        	
    }
    
    /**
     * check if the atom is saturated
     * @param ac
     * @param atom
     * @return
     */
    private static boolean isSaturate(IAtomContainer ac, IAtom atom) {
    	double bondOrderSum = ac.getBondOrderSum(atom);
		int hcount = atom.getImplicitHydrogenCount() == CDKConstants.UNSET ? 0 : atom.getImplicitHydrogenCount();
		int charge = atom.getFormalCharge() == CDKConstants.UNSET ? 0 : atom.getFormalCharge();
		int radical = atom.getProperty(additionalConstants.RADICAL) == null ? 0 : atom.getProperty(additionalConstants.RADICAL);

		
		// calculate the valence
		double cVal = bondOrderSum - charge + hcount + radical;
    	return ElementCalculation.calculateValence(Element.valueOfIgnoreCase(atom.getSymbol())) == cVal;
    	
    }
    
    
    static int initial(final IAtomContainer ac, Matching2 m, final BitSet s) {

        int nMatched = 0;

        for (int v = s.nextSetBit(0); v >= 0; v = s.nextSetBit(v + 1)) {
        	IAtom atom = ac.getAtom(v);
            // skip if already matched
            if (m.matched(v))
                continue;

            // find a single edge which is not matched and match it
            for (IBond e : ac.getConnectedBondsList(atom)) {
                int w = e.getOther(atom).getIndex();
                if ((!e.getOrder().equals(IBond.Order.SINGLE)) && m.unmatched(w) && s.get(w)) {
                    m.match(v, w);
                    nMatched += 2;
                    break;
                }
            }
        }

        return nMatched;
    }
    
    /**
     * Determine the set of atoms that are available to have a double-bond.
     *
     * @param graph adjacent list representation
     * @param atoms array of atoms
     * @param bonds map of atom indices to bonds
     * @return atoms that can require a double-bond
     */
    private static BitSet available(int[][] graph, IAtom[] atoms, EdgeToBondMap bonds) {

        final BitSet available = new BitSet();

        // for all atoms, select those that require a double-bond
        ATOMS: for (int i = 0; i < atoms.length; i++) {

            final IAtom atom = atoms[i];

            // preconditions
            if (atom.getAtomicNumber() == null)
                throw new IllegalArgumentException("atom " + (i + 1) + " had unset atomic number");
            if (atom.getFormalCharge() == null)
                throw new IllegalArgumentException("atom " + (i + 1) + " had unset formal charge");
            if (atom.getImplicitHydrogenCount() == null)
                throw new IllegalArgumentException("atom " + (i + 1) + " had unset implicit hydrogen count");

            if (!atom.isAromatic()) continue;

            // count preexisting pi-bonds, a higher bond order causes a skip
            int nPiBonds = 0;
            for (final int w : graph[i]) {
                IBond.Order order = bonds.get(i, w).getOrder();
                if (order == DOUBLE) {
                    nPiBonds++;
                } else if (order.numeric() > 2) {
                    continue ATOMS;
                }
            }	

            // check if a pi bond can be assigned
            final int element = atom.getAtomicNumber();
            final int charge = atom.getFormalCharge();
            final int valence;
            //a piBons has been already attributed, so the atom is not available
            if (nPiBonds > 0)
            	continue;
            if (atom.getProperty("bondOrderSum2") == null)
            	valence = graph[i].length + atom.getImplicitHydrogenCount() + nPiBonds;
            else{
            	valence = (int)atom.getProperty("bondOrderSum2") + atom.getImplicitHydrogenCount();
            }
            
            if (available(element, charge, valence)) {
                available.set(i);
            }
        }
        return available;
    }

    
    /**
     * Determine if the specified element with the provided charge and valance
     * requires a pi bond?
     *
     * @param element atomic number >= 0
     * @param charge  formal charge
     * @param valence bonded electrons
     * @return a double-bond is required
     */
    private static boolean available(final int element, final int charge, final int valence) {
        // higher atomic number elements aren't likely to be found but
        // we have them for rare corner cases (tellurium).
        // Germanium, Silicon, Tin and Antimony are a bit bonkers...
        switch (Elements.ofNumber(element)) {
            case Boron:
                if (charge == 0 && valence <= 2) return true;
                if (charge == -1 && valence <= 3) return true;
                break;
            case Carbon:
            case Silicon:
            case Germanium:
            case Tin:
                if (charge == 0 && valence <= 3) return true;
                break;
            case Nitrogen:
            case Phosphorus:
            case Arsenic:
            case Antimony:
                if (charge == 0) return valence <= 2 || valence == 4;
                if (charge == 1) return valence <= 3;
                break;
            case Oxygen:
            case Sulfur:
            case Selenium:
            case Tellurium:
                // valence of three or five are really only for sulphur but
                // are applied generally to all of group eight for simplicity
                if (charge == 0) return valence <= 1 || valence == 3 || valence == 5;
                if (charge == 1) return valence <= 2 || valence == 4;
                break;
        }

        return false;
    }
    
    
    static BitSet buildSet(IAtomContainer ac, BitSet aromatic) {

        BitSet undecided = new BitSet(ac.getAtomCount());

        for (int v = 0; v < ac.getAtomCount(); v++) {
            IAtom atom = ac.getAtom(v);
        	if (atom.isAromatic()) {
                aromatic.set(v);
                if (!predetermined2(ac, v))
                    undecided.set(v);
            }
        }

        return undecided;
    }
    
    static boolean predetermined(IAtomContainer ac, int v) {

        IAtom a = ac.getAtom(v);
        int bondOrderSum = (int) ac.getBondOrderSum(a);

        int q = a.getFormalCharge();
        int deg = bondOrderSum + a.getImplicitHydrogenCount();
        //this valence exclude the number of implicit Hydrogen
        int valence = bondOrderSum - q;

        if (valence > bondOrderSum) {
        	for (IBond b : ac.getConnectedBondsList(a)) {
        		if (b.getOrder() == IBond.Order.DOUBLE) {
                    if (q == 0 && (a.getSymbol().equals("N") || (a.getSymbol().equals("S") && deg > 3))
                            )
                        return false;
                    return true;
                }
                // triple or quadruple bond - we don't need to assign anymore p electrons
                else if (b.getOrder().numeric() > 2) {
                    return true;
                }
        	}
        }
        
     // no pi bonds does the degree and charge indicate that
        // there can be no other pi bonds
        switch (a.getSymbol()) {
            case "B":
                return (q == 0) && deg == 3;
            case "C":
                return (q == 1 || q == -1) && deg == 3;
            case "Si":
            case "Ge":
                return q < 0;
            case "N":
            case "P":
            case "As":
            case "Sb":
                if (q == 0)
                    return deg == 3 || deg > 4;
                else if (q == 1)
                    return deg > 3;
                else
                    return true;
            case "O":
            case "S":
            case "Se":
            case "Te":
                if (q == 0)
                    return deg == 2 || deg == 4 || deg > 5;
                else if (q == -1 || q == +1)
                    return deg == 3 || deg == 5 || deg > 6;
                else
                    return false;
        }

        return false;
    }
    
    static boolean predetermined2(IAtomContainer ac, int v) {

        IAtom a = ac.getAtom(v);
        int bondOrderSum = (int) ac.getBondOrderSum(a);

        int q = a.getFormalCharge();
        int deg = bondOrderSum + a.getImplicitHydrogenCount();
        //this valence exclude the number of implicit Hydrogen
        int valence = bondOrderSum - q;

        if (valence > bondOrderSum) {
        	for (IBond b : ac.getConnectedBondsList(a)) {
        		if (b.getOrder() == IBond.Order.DOUBLE) {
                    if (q == 0 && (a.getSymbol().equals("N") || (a.getSymbol().equals("S") && deg > 3))
                            )
                        return false;
                    return true;
                }
                // triple or quadruple bond - we don't need to assign anymore p electrons
                else if (b.getOrder().numeric() > 2) {
                    return true;
                }
        	}
        }
        
     // no pi bonds does the degree and charge indicate that
        // there can be no other pi bonds
        switch (a.getSymbol()) {
            case "B":
                return (q == 0) && deg == 3;
            case "C":
                return (q == 1 || q == -1) && deg == 3;
            case "Si":
            case "Ge":
                return q < 0;
            case "N":
            	return deg == 3;
            case "P":
            	return deg == 3;
            case "As":
            case "Sb":
                if (q == 0)
                    return deg == 3 || deg > 4;
                else if (q == 1)
                    return deg > 3;
                else
                    return true;
            case "O":
            	return deg == 2;
            case "S":
            	return deg == 2;
            case "Se":
            case "Te":
                if (q == 0)
                    return deg == 2 || deg == 4 || deg > 5;
                else if (q == -1 || q == +1)
                    return deg == 3 || deg == 5 || deg > 6;
                else
                    return false;
        }

        return false;
    }   
    /**
     * UNUSED
     * @param ac
     */
    /*
    private static void pseudoRandomlyDoubleBondAttribution(IAtomContainer ac) {
    	CycleFinder cf = Cycles.mcb();
    	try {
    		Cycles   cycles = cf.find(ac);
    		IRingSet rings  = cycles.toRingSet();
    		List<IRingSet> systems = RingPartitioner.partitionRings(rings);
    		for (IRingSet system : systems) {
    			boolean isAromatic = false;
    			int bestCount = 999999;
    			List<IRing> bestRings = new ArrayList<IRing>();
    			Set<IAtom> atoms = new HashSet<IAtom>();
    			//get problematic ring system (odd atom number)
    			for (IAtomContainer ring : system.atomContainers()) {
    				int count = system.getConnectedRings((IRing) ring).getAtomContainerCount();
    				if (count < bestCount) {
    					bestCount = count;
    					bestRings.clear();
    					bestRings.add((IRing) ring);

    				}
    				else if (count == bestCount) {
    					bestRings.add((IRing) ring);
    				}
    				for (IAtom atom : ring.atoms()) {
    					atoms.add(atom);
    					if (atom.isAromatic())
    						isAromatic = true;
    				}
    			}
    			if (isAromatic && (atoms.size() & 0x1) == 0x1) {
    				IRing best = null;
    				int n = 99999;
    				boolean valid = false;
    				for (IRing ring : bestRings) {
    					int cpt = ring.getAtomCount();
    					if (cpt < n) {
    						for (IBond b : ring.bonds()) {
    							if (b.getBegin().getSymbol().equals("C") && b.getEnd().getSymbol().equals("C") &&
    									ac.getBondOrderSum(b.getBegin()) == 2 && ac.getBondOrderSum(b.getEnd()) == 2) {
    								valid = true;
    								break;
    							}
    						}
    						if (valid) {
    							n = cpt;
    							best = ring;
    						}
    					}
    				}
    				for (IBond b : best.bonds()) {
    					if (b.getBegin().getSymbol().equals("C") && b.getEnd().getSymbol().equals("C") &&
    							ac.getBondOrderSum(b.getBegin()) == 2 && ac.getBondOrderSum(b.getEnd()) == 2) {
    						b.setOrder(IBond.Order.DOUBLE);
    						break;
    					}
    				}
    			}
    		} 
    	}catch (Intractable e) {
    		// ignore error - MCB should never be intractable
    	}
    
    

 
    private static boolean forceKekulize(IAtomContainer ac) throws CDKException {
    	BitSet aromatic = new BitSet();
        BitSet subset = buildSet(ac, aromatic);
    	if (hasOddCardinality(subset)){
    		//throw new CDKException("a valid kekulé structure could not be assigned");
    		return false;
    	}
    	final Matching2 m = Matching2.empty(ac.getAtomCount());
    	int nMatched = initial(ac, m, subset);
    	assign(ac, subset, aromatic, m);
    	return true;
    }
    
    
    private static boolean hasOddCardinality(BitSet s) {
        return (s.cardinality() & 0x1) == 1;
    }
    
    // invariant, m is a perfect matching
    private static IAtomContainer assign(IAtomContainer ac, BitSet subset, BitSet aromatic, Matching2 m) {
         for (int u = aromatic.nextSetBit(0); u >= 0; u = aromatic.nextSetBit(u + 1)) {
        	 IAtom atom = ac.getAtom(u);
             for (IBond e : ac.getConnectedBondsList(atom)) {
                 int v = e.getOther(atom).getIndex();
                 if (v < u) {
                	 if (e.isAromatic()) {
                		 if (e.getOrder().equals(IBond.Order.SINGLE)) {
                    		 if (subset.get(u) && m.other(u) == v) {
                                 e.setOrder(IBond.Order.DOUBLE);
                             }
                             else if (aromatic.get(v)) {
                            	 e.setOrder(IBond.Order.UNSET);
                             }
                             else {
                            	 e.setOrder(IBond.Order.SINGLE);
                             }
                    	 }
                		 else {
                			 if (subset.get(u) && m.other(u) == v) {
                            	 e.setOrder(IBond.Order.DOUBLE);
                             }
                             else if (aromatic.get(u) && aromatic.get(v)) {
                            	 e.setOrder(IBond.Order.UNSET);
                             }
                		 }
                	 }
                 }
             }
         }
         return ac;
    */
}

