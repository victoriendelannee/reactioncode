package com.nih.fragments;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.ringsearch.RingSearch;

import com.nih.util.tools;

public class Murcko {
	private static final String IS_SIDECHAIN_ATOM    = "sidechainatom";
    private static final String IS_LINKER_ATOM       = "linkeratom";
    private static final String IS_CONNECTED_TO_RING = "rcon";
    private static final String IS_IN_MURCKO_SCAFFOLD    = "murcko";
    
    static Map<IAtom,List<IAtom>> connections; 
    
    public static void initMarkAtoms(IAtomContainer ac) throws CDKException {
        // identify rings
		//Method 1
    	//RingSearch ringSearch = new RingSearch(ac);
    	//IAtomContainer r = ringSearch.ringFragments();
    	
    	//for (IBond bond : r.bonds())
            //bond.setFlag(CDKConstants.ISINRING, true);
    	
    	//Method 2
    	//Cycles allCycles = Cycles.sssr(ac);
		//IRingSet r = allCycles.toRingSet();
    	
    	//Method 3
    	AllRingsFinder arf = new AllRingsFinder().usingThreshold(AllRingsFinder.Threshold.None);
    	IRingSet r = arf.findAllRings(ac);
    	
//        // manually flag ring bonds
        for (IAtomContainer ar : r.atomContainers()) {
            for (IBond bond : ar.bonds())
                bond.setFlag(CDKConstants.ISINRING, true);
        }

        int cpt = 0;
        for (IAtom atom : ac.atoms()) {
            atom.setProperty(IS_LINKER_ATOM, false);
            atom.setProperty(IS_SIDECHAIN_ATOM, false);
            atom.setProperty(IS_CONNECTED_TO_RING, false);
            atom.setProperty("uID", cpt);
            cpt++;
        }
    }
    

    
    public static void markLinkers(IAtomContainer atomContainer) {
    	//to manage case where a non ring bond link 2 ring system

    	// first we check for single atoms between rings - these are linker atoms
        // this is also the place where we need to check for something like PhC(C)Ph
        // since the central atom is a single atom between rings, but also has a non
        // ring attachment
        for (IAtom atom : atomContainer.atoms()) {
        	//case where a non ring bond link 2 ring system
        	List<IBond> conBonds = atomContainer.getConnectedBondsList(atom);
        	if (atom.getFlag(CDKConstants.ISINRING)) {
        		for (IBond conBond : conBonds) {
                    if (!conBond.getFlag(CDKConstants.ISINRING)) {
                        if (conBond.getBegin().getFlag(CDKConstants.ISINRING) && conBond.getEnd().getFlag(CDKConstants.ISINRING))
                        	atom.setProperty(IS_LINKER_ATOM, true);
                    }
                }
            }
            List<IAtom> conatoms = atomContainer.getConnectedAtomsList(atom);
            if (conatoms.size() == 1) continue; // this is actually a terminal atom and so is a side chain
            int nRingAtom = 0;
            for (IAtom conatom : conatoms) {
                if (conatom.getFlag(CDKConstants.ISINRING)) {
                    nRingAtom++;
                }
            }
            if (nRingAtom > 0) atom.setProperty(IS_CONNECTED_TO_RING, true);
            if (nRingAtom >= 2) atom.setProperty(IS_LINKER_ATOM, true);
        }

        // now lets look at linker paths
        for (IAtom atom1 : atomContainer.atoms()) {
            if (atom1.getFlag(CDKConstants.ISINRING) || !(Boolean) atom1.getProperty(IS_CONNECTED_TO_RING)) continue;
            for (IAtom atom2 : atomContainer.atoms()) {
                if (atom2.getFlag(CDKConstants.ISINRING) || !(Boolean) atom2.getProperty(IS_CONNECTED_TO_RING))
                    continue;

                if (atom1.equals(atom2)) continue;

                // ok, get paths between these two non-ring atoms. Each of these atoms
                // is connected to a ring atom, and so if the atoms between these atoms
                // not ring atoms, this is a linker path
                List<List<IAtom>> paths = PathTools.getAllPaths(atomContainer, atom1, atom2);

                for (List<IAtom> path : paths) {
                    boolean allNonRing = true;
                    for (IAtom atom : path) {
                        if (atom.getFlag(CDKConstants.ISINRING)) {
                            allNonRing = false;
                            break;
                        }
                    }
                    if (allNonRing) { // mark them as linkers
                        for (IAtom atom : path)
                            atom.setProperty(IS_LINKER_ATOM, true);
                    }
                }
            }
        }
    }

	//CDK implementation: a non ring bond connected 2 ring System in this present code is not considered as a linker
    public static void markLinkers2(IAtomContainer atomContainer) {
        // first we check for single atoms between rings - these are linker atoms
        // this is also the place where we need to check for something like PhC(C)Ph
        // since the central atom is a single atom between rings, but also has a non
        // ring attachment
        for (IAtom atom : atomContainer.atoms()) {
            if (atom.getFlag(CDKConstants.ISINRING)) continue; // only need to look at non-ring atoms
            List<IAtom> conatoms = atomContainer.getConnectedAtomsList(atom);
            if (conatoms.size() == 1) continue; // this is actually a terminal atom and so is a side chain
            int nRingAtom = 0;
            for (IAtom conatom : conatoms) {
                if (conatom.getFlag(CDKConstants.ISINRING)) {
                    nRingAtom++;
                }
            }
            if (nRingAtom > 0) atom.setProperty(IS_CONNECTED_TO_RING, true);
            if (nRingAtom >= 2) atom.setProperty(IS_LINKER_ATOM, true);
        }

        // now lets look at linker paths
        for (IAtom atom1 : atomContainer.atoms()) {
            if (atom1.getFlag(CDKConstants.ISINRING) || !(Boolean) atom1.getProperty(IS_CONNECTED_TO_RING)) continue;
            for (IAtom atom2 : atomContainer.atoms()) {
                if (atom2.getFlag(CDKConstants.ISINRING) || !(Boolean) atom2.getProperty(IS_CONNECTED_TO_RING))
                    continue;

                if (atom1.equals(atom2)) continue;

                // ok, get paths between these two non-ring atoms. Each of these atoms
                // is connected to a ring atom, and so if the atoms between these atoms
                // not ring atoms, this is a linker path
                List<List<IAtom>> paths = PathTools.getAllPaths(atomContainer, atom1, atom2);

                for (List<IAtom> path : paths) {
                    boolean allNonRing = true;
                    for (IAtom atom : path) {
                        if (atom.getFlag(CDKConstants.ISINRING)) {
                            allNonRing = false;
                            break;
                        }
                    }
                    if (allNonRing) { // mark them as linkers
                        for (IAtom atom : path)
                            atom.setProperty(IS_LINKER_ATOM, true);
                    }
                }
            }
        }
    }

    public static void markSideChains(IAtomContainer atomContainer) {
        for (IAtom atom : atomContainer.atoms()) {
            if (!isring(atom) && !islinker(atom)) atom.setProperty(IS_SIDECHAIN_ATOM, true);
        }
    }
    
    public static void markMurckoScaffold(IAtomContainer atomContainer) {
        for (IAtom atom : atomContainer.atoms()) {
            if (isring(atom) || islinker(atom)) atom.setProperty(IS_IN_MURCKO_SCAFFOLD, true);
            else atom.setProperty(IS_IN_MURCKO_SCAFFOLD, false);
        }
    }
    
    /**
     * Computes the Murcko Scaffold for the provided molecule in linear time.
     * Note the return value contains the same atoms/bonds as in the input
     * and an additional clone and valence adjustments may be required.
     *
     * @param mol the molecule
     * @return the atoms and bonds in the scaffold
     */
    public static IAtomContainer scaffold(final IAtomContainer mol) {

    	connections = new HashMap<IAtom,List<IAtom>>();
    	
        // Old AtomContainer IMPL, cannot work with this
        if (!mol.isEmpty() && mol.getAtom(0).getContainer() == null)
            return null;

        Deque<IAtom> queue = new ArrayDeque<>();
        int[] bcount = new int[mol.getAtomCount()];

        // Step 1. Mark and queue all terminal (degree 1) atoms
        for (IAtom atom : mol.atoms()) {
            int numBonds = atom.getBondCount();
            bcount[atom.getIndex()] = numBonds;
            if (numBonds == 1)
                queue.add(atom);
        }

        // Step 2. Iteratively remove terminal atoms queuing new atoms
        //         as they become terminal
        while (!queue.isEmpty()) {
            IAtom atom = queue.poll();
            if (atom == null)
                continue;
            bcount[atom.getIndex()] = 0;
            for (IBond bond : atom.bonds()) {
                IAtom nbr = bond.getOther(atom);
                bcount[nbr.getIndex()]--;
                if (bcount[nbr.getIndex()] == 1)
                    queue.add(nbr);
                else {
                	if (((Boolean) atom.getProperty(IS_LINKER_ATOM) == true || (Boolean) atom.getProperty(IS_SIDECHAIN_ATOM) == true) && 
                			nbr.getFlag(CDKConstants.ISINRING)) {
                		if (!connections.containsKey(nbr)) {
                			List<IAtom> l = new ArrayList<IAtom>();
                			l.add(atom);
                			connections.put(nbr, l);
                		}
                		else {
                			List<IAtom> l = connections.get(nbr);
                			l.add(atom);
                			connections.put(nbr, l);
                		}
                	}
                }
            }
        }

        // Step 3. Copy out the atoms/bonds that are part of the Murcko
        //         scaffold
        IAtomContainer scaffold = mol.getBuilder().newAtomContainer();
        for (int i = 0; i < mol.getAtomCount(); i++) {
            IAtom atom = mol.getAtom(i);
            if (bcount[i] > 0)
                scaffold.addAtom(atom);
        }
        for (int i = 0; i < mol.getBondCount(); i++) {
            IBond bond = mol.getBond(i);
            if (bcount[bond.getBegin().getIndex()] > 0 &&
                bcount[bond.getEnd().getIndex()] > 0)
                scaffold.addBond(bond);
        }
        return scaffold;
    }
    
    
    private static boolean isring(IAtom atom) {
        return atom.getFlag(CDKConstants.ISINRING);
    }

    private static boolean islinker(IAtom atom) {
        Boolean o = (Boolean) atom.getProperty(IS_LINKER_ATOM);
        return o != null && o;
    }

    private boolean issidechain(IAtom atom) {
        Boolean o = (Boolean) atom.getProperty(IS_SIDECHAIN_ATOM);
        return o != null && o;
    }

    private boolean islinker(IBond bond) {
        return islinker(bond.getBegin()) || islinker(bond.getEnd());
    }
    
    private boolean isZeroAtomLinker(IBond bond) {
        boolean isRingBond = bond.getFlag(CDKConstants.ISINRING);
        return isring(bond.getBegin()) && isring(bond.getEnd()) && !isRingBond;
    }
    
    private boolean hasframework(IAtomContainer atomContainer) {
        boolean hasLinker = false;
        boolean hasRing = false;
        for (IAtom atom : atomContainer.atoms()) {
            if (islinker(atom)) hasLinker = true;
            if (isring(atom)) hasRing = true;
            if (hasLinker && hasRing) break;
        }

        // but two rings may be connected by a single bond
        // in which case, the atoms of the bond are not
        // linker atoms, but the bond itself is a (pseudo) linker bond
        for (IBond bond : atomContainer.bonds()) {
            if (isZeroAtomLinker(bond)) {
                hasLinker = true;
                break;
            }
        }
        return hasLinker && hasRing;
    }



	public static Map<IAtom, List<IAtom>> getConnections() {
		return connections;
	}
    
}
