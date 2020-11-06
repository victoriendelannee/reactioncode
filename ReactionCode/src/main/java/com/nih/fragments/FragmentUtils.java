package com.nih.fragments;
/* Copyright (C) 2010  Rajarshi Guha <rajarshi.guha@gmail.com>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */

import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Helper methods for fragmentation algorithms.
 * 
 * Most of these methods are specific to the fragmentation algorithms
 * in this package and so are protected. In general, these methods will
 * not be used by the rest of the API or by other users of the library.
 *
 * @author Rajarshi Guha
 * @cdk.module fragment
 */
public class FragmentUtils {

	/**
	 * Non destructively split a molecule into two parts at the specified bond.
	 *
	 * Note that if a ring bond is specified, the resultant list will contain
	 * the opened ring twice.
	 *
	 * @param atomContainer The molecule to split
	 * @param bond The bond to split at
	 * @return A list containing the two parts of the molecule
	 */
	public static List<IAtomContainer> splitMolecule(IAtomContainer atomContainer, IBond bond) {
		List<IAtomContainer> ret = new ArrayList<IAtomContainer>();

		for (IAtom atom : bond.atoms()) {

			// later on we'll want to make sure that the fragment doesn't contain
			// the bond joining the current atom and the atom that is on the other side
			IAtom excludedAtom = null;
			if (atom.equals(bond.getAtom(0)))
				excludedAtom = bond.getAtom(1);
			else
				excludedAtom = bond.getAtom(0);

			List<IBond> part = new ArrayList<IBond>();
			part.add(bond);
			part = traverse(atomContainer, atom, part);

			// at this point we have a partition which contains the bond we
			// split. This partition should actually 2 partitions:
			// - one with the splitting bond
			// - one without the splitting bond
			// note that this will lead to repeated fragments when we  do this
			// with adjacent bonds, so when we gather all the fragments we need
			// to check for repeats
			IAtomContainer partContainer;
			partContainer = makeAtomContainer(atom, part, excludedAtom);

			// by checking for more than 2 atoms, we exclude single bond fragments
			// also if a fragment has the same number of atoms as the parent molecule,
			// it is the parent molecule, so we exclude it.
			if (partContainer.getAtomCount() > 0 && partContainer.getAtomCount() != atomContainer.getAtomCount())
				ret.add(partContainer);

			//part.remove(0);
			//partContainer = makeAtomContainer(atom, part, excludedAtom);
			//if (partContainer.getAtomCount() > 2 && partContainer.getAtomCount() != atomContainer.getAtomCount())
			//	ret.add(partContainer);
		}

		return ret;
	}
	
	/**
	 * Non destructively split a molecule into two parts at the specified bond and return desire part.
	 *
	 * Note that if a ring bond is specified, the resultant list will contain
	 * the opened ring twice.
	 *
	 * @param atomContainer The molecule to split
	 * @param bond The bond to split at
	 * @return A list containing the two parts of the molecule
	 */
	public static IAtomContainer splitMolecule(IAtomContainer atomContainer, IBond bond, IAtom excludedAtom, boolean addSplittingBond) {
			IAtom other = bond.getOther(excludedAtom);

			IAtomContainer partContainer = other.getBuilder().newInstance(IAtomContainer.class);
//			partContainer.setAtom(0, other);
			partContainer.addAtom(other);
			
			
			Set<IAtom> sphere = new HashSet<IAtom>();
			sphere.add(other);
			//set as visited to exclude the desire part of the molecule
			bond.setFlag(CDKConstants.VISITED, true);
			makeAtomContainer(sphere, atomContainer, partContainer);

			if (addSplittingBond == true) {
				Bond b = new Bond();
				b.setAtom(other, 0);
				b.setAtom(new Atom("*"), 1);
				b.setOrder(bond.getOrder());
				b.setStereo(bond.getStereo());
				partContainer.addAtom(b.getAtom(1));
//				partContainer.setAtom(partContainer.getAtomCount()-1, b.getAtom(1));
				partContainer.addBond(b);
			}

		return partContainer;
	}
	
	/**
	 * Second method (slowest)
	 * Non destructively split a molecule into two parts at the specified bond and return desire part.
	 *
	 * Note that if a ring bond is specified, the resultant list will contain
	 * the opened ring twice.
	 *
	 * @param atomContainer The molecule to split
	 * @param bond The bond to split at
	 * @return A list containing the two parts of the molecule
	 */
	public static IAtomContainer splitMolecule2(IAtomContainer atomContainer, IBond bond, IAtom excludedAtom, boolean addSplittingBond) {
			List<IBond> part = new ArrayList<IBond>();
			part.add(bond);
			IAtom other = bond.getOther(excludedAtom);
			part = traverse(atomContainer, other, part);

			// at this point we have a partition which contains the bond we
			// split. This partition should actually 2 partitions:
			// - one with the splitting bond
			// - one without the splitting bond
			// note that this will lead to repeated fragments when we  do this
			// with adjacent bonds, so when we gather all the fragments we need
			// to check for repeats
			
			IAtomContainer partContainer;
			partContainer = makeAtomContainer(other, part, excludedAtom);

			if (addSplittingBond == true) {
				Bond b = new Bond();
				b.setAtom(other, 0);
				b.setAtom(new Atom("*"), 1);
				b.setOrder(bond.getOrder());
				b.setStereo(b.getStereo());
				partContainer.addAtom(b.getAtom(1));
//				partContainer.setAtom(partContainer.getAtomCount()-1, b.getAtom(1));
				partContainer.addBond(b);
			}

		return partContainer;
	}

	public static List<IBond> traverse(IAtomContainer atomContainer, IAtom atom, List<IBond> bondList) {
		List<IBond> connectedBonds = atomContainer.getConnectedBondsList(atom);
		for (IBond aBond : connectedBonds) {
			if (bondList.contains(aBond)) continue;
			bondList.add(aBond);
			IAtom[] atoms = getIAtomArray(aBond);
			IAtom nextAtom = null;
			if (atoms[0].equals(atom))
	            nextAtom = atoms[1];
	        else if (atoms[1].equals(atom))
	            nextAtom = atoms[0];
			if (atomContainer.getConnectedBondsCount(nextAtom) == 1) continue;
			traverse(atomContainer, nextAtom, bondList);
		}
		return bondList;
	}
	
	protected static IAtom[] getIAtomArray(IBond bond) {
		IAtom[] atoms = new IAtom[2];
        atoms[0] = bond.getAtom(0);
        atoms[1] = bond.getAtom(1);
        return atoms;
	}
	
	protected static IAtomContainer makeAtomContainer(IAtom atom, List<IBond> parts, IAtom excludedAtom) {
		IAtomContainer partContainer = atom.getBuilder().newInstance(IAtomContainer.class);
		partContainer.addAtom(atom);
		for (IBond aBond : parts) {
			for (IAtom bondedAtom : aBond.atoms()) {
				if (!bondedAtom.equals(excludedAtom) && !partContainer.contains(bondedAtom))
					partContainer.addAtom(bondedAtom);
			}
			if (!aBond.contains(excludedAtom)) partContainer.addBond(aBond);
		}
		return partContainer;
	}
	
	public static List<IBond> getAllConnectedBonds(Set<IAtom> sphere, IAtomContainer atomContainer, List<IBond> bondList, Set<IAtom> visited) {
		IAtom nextAtom;
		Set<IAtom> newSphere = new HashSet<IAtom>();

		for (IAtom atom : sphere) {
			List<IBond> connectedBonds = atomContainer.getConnectedBondsList(atom);
			for (IBond bond : connectedBonds) {
				nextAtom = bond.getConnectedAtom(atom);
				if (!bondList.contains(bond)) {
					bondList.add(bond);
				}
				if (!visited.contains(nextAtom)) {
					if (!sphere.contains(nextAtom)) newSphere.add(nextAtom);
					visited.add(nextAtom);
				}
			}
		}
		if (newSphere.size() > 0) {
			getAllConnectedBonds(newSphere, atomContainer, bondList, visited);
		}
		return bondList;
	}
	
	public static IAtomContainer makeAtomContainer(IAtom atom, List<IBond> parts) {
		IAtomContainer partContainer = atom.getBuilder().newInstance(IAtomContainer.class);
		partContainer.addAtom(atom);
		for (IBond aBond : parts) {
			for (IAtom bondedAtom : aBond.atoms()) {
				if (!partContainer.contains(bondedAtom))
					partContainer.addAtom(bondedAtom);
			}
			partContainer.addBond(aBond);
		}
		return partContainer;
	}
	
	private static void makeAtomContainer(Set<IAtom> sphere, IAtomContainer atomContainer, IAtomContainer partContainer) {
		IAtom nextAtom;
		Set<IAtom> newSphere = new HashSet<IAtom>();

		for (IAtom atom : sphere) {
			List<IBond> bonds = atomContainer.getConnectedBondsList(atom);
			for (IBond bond : bonds) {
				nextAtom = bond.getConnectedAtom(atom);
				if (!nextAtom.getFlag(CDKConstants.VISITED)) {
					if (!sphere.contains(nextAtom)) newSphere.add(nextAtom);
					partContainer.addAtom(nextAtom);
					nextAtom.setFlag(CDKConstants.VISITED, true);
				}
				if (!partContainer.contains(bond)) 
					partContainer.addBond(bond);
			}
		}
		if (newSphere.size() > 0) {
			makeAtomContainer(newSphere, atomContainer, partContainer);
		}
	}

	public static IAtomContainerSet makeAtomContainerSet(IAtomContainer ac) {
		IAtomContainerSet set = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		
		Map<IChemObject, IStereoElement> seMap = new HashMap<IChemObject, IStereoElement>();
		
		for (IStereoElement se : ac.stereoElements()) {
			seMap.put(se.getFocus(), se);
		}
		
		for (IAtom a : ac.atoms()) {
			if (!a.getFlag(CDKConstants.VISITED)) {
				IAtomContainer partContainer;
				if (ac instanceof QueryAtomContainer)
					partContainer = new QueryAtomContainer(DefaultChemObjectBuilder.getInstance());
				else
					partContainer = a.getBuilder().newInstance(IAtomContainer.class);
				partContainer.addAtom(a);
				a.setFlag(CDKConstants.VISITED, true);

				Set<IAtom> sphere = new HashSet<IAtom>();
				sphere.add(a);
				makeAtomContainer(sphere, ac, partContainer);
				
				if (!seMap.isEmpty())
					partContainer.setStereoElements(getStereoElement(partContainer, seMap));
				set.addAtomContainer(partContainer);
				
			}
		}
		return set;
	}
	
	private static List<IStereoElement> getStereoElement(IAtomContainer ac, Map<IChemObject, IStereoElement> seMap) {
		List<IStereoElement> result = new ArrayList<IStereoElement>();
		
		for (IAtom a : ac.atoms()) {
			if (seMap.containsKey(a))
				result.add(seMap.get(a));
		}
		
		for (IBond b : ac.bonds()) {
			if (seMap.containsKey(b))
				result.add(seMap.get(b));
		}
		return result;
	}
	
	public static void resetVisitedFlags(IAtomContainerSet set) {
		for (IAtomContainer ac : set.atomContainers()) {
			resetVisitedFlags(ac);
		}
	}
	
	public static void resetVisitedFlags(IAtomContainer ac) {
		for (IAtom a : ac.atoms()) {
			a.setFlag(CDKConstants.VISITED, false);
		}
		
		for (IBond b : ac.bonds()) {
			b.setFlag(CDKConstants.VISITED, false);
		}
	}
}
