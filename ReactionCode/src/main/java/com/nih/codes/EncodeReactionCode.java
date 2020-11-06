/**
 * 
 */
package com.nih.codes;

import static org.openscience.cdk.CDKConstants.ATOM_ATOM_MAPPING;
import static com.nih.reaction.additionalConstants.BOND_CHANGE_INFORMATION;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.Map.Entry;


import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.Bond;
import org.openscience.cdk.tools.periodictable.PeriodicTable;

import org.openscience.cdk.smsd.helper.BondEnergy;

import com.nih.reaction.additionalConstants;
import com.nih.tools.ElementCalculation;


/**
 * @author delanneev
 *
 */
public class EncodeReactionCode {

	private IAtomContainer mol;
	private HashMap<Integer, String> connectionTableAlphabet = new HashMap<Integer, String>();
	private String reactionCode = null;
	
	//parameters
	private boolean hybridization = false;
	private boolean stereochemistry = true;
	private boolean charge = true;
	private boolean bondType = true;
	private boolean repetition = true;
	
	private Map<Integer, Set<IAtom>> codeLayers; 
	private List<IAtom> conflictedAtoms;
	
	private Set<Integer> atomInReactionCenterNotInProducts;

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}
	
	public EncodeReactionCode() {
		atomInReactionCenterNotInProducts = new HashSet<Integer>();
	}
	
	public EncodeReactionCode(Set<Integer> atomInReactionCenterNotInProducts) {
		this.atomInReactionCenterNotInProducts = atomInReactionCenterNotInProducts;
		//make table of correspondence to connect atom for reaction code generation
//		connectionTableAlphabet = makeConnectionTableAlphabet();
	}

	/**
	 * 
	 * @param atomInRC the bonds in the reaction center
	 * @param reactants to get the hybridization because the pseudoMolecule contains the hybridization of the products
	 * @param pseudoMolecule
	 * @param numberOfRepetitions the number of time each atom is duplicate ex: c1(c2ccc(Br)cc2)ccc(Br)cc1.CC1=CC(=CC=C1)NC2=CC=CC=C2>>N(c1ccc(c2ccc(N(c3cccc(C)c3)c3ccccc3)cc2)cc1)(c1cccc(C)c1)c1ccccc1
	 * @return
	 */
	public Map<String,String> generateReactionCode(Set<IAtom> atomInRC, IAtomContainerSet reactants,
			IAtomContainer pseudoMolecule, Map<String, Integer> numberOfRepetitions) {

		mol = pseudoMolecule;
		reactionCode = null;
		
		//make table of correspondence to connect atom for reaction code generation
		connectionTableAlphabet = makeConnectionTableAlphabet();
		
		//generate index
		Map<String, IAtom> indexAtomReactants = indexAtomIAtomContainerSet(reactants);
		Map<String, IBond> indexBondReactants = indexAtomIBondContainerSet(reactants);
		reinitializedVisitedConstants(pseudoMolecule);

		//generate all individual atom and bonds code
		codeLayers = new LinkedHashMap<Integer, Set<IAtom>>();
		codeLayers.put(0, atomInRC);
		generateAllAtomAndBondCodes(atomInRC, indexAtomReactants, indexBondReactants, pseudoMolecule, codeLayers, 
				new ArrayList<IAtom>(), 1, pseudoMolecule.getAtomCount());
		
		return generate(pseudoMolecule, codeLayers, numberOfRepetitions);
	}



	/**
	 * @param pseudoMol
	 * @param codeLayers
	 * @param numberOfRepetitions
	 * @return
	 */
	private Map<String,String> generate(IAtomContainer pseudoMol, Map<Integer, Set<IAtom>> codeLayers,
			Map<String, Integer> numberOfRepetitions) {
		int position = 1;
		for (Entry<Integer, Set<IAtom>> e : codeLayers.entrySet()) {
			int depth = e.getKey();
			List<IAtom> sortedAtomsByAtomCode = new ArrayList<IAtom>(e.getValue());
			List<IAtom> sortedAtomsByABCode = new ArrayList<IAtom>(e.getValue());
			Collections.sort(sortedAtomsByAtomCode, new CompareByAtomCode());
			Collections.sort(sortedAtomsByABCode, new CompareByABCode());
			
			for (int i = 0; i < sortedAtomsByAtomCode.size(); i++) {
				conflictedAtoms = new ArrayList<IAtom>();
				IAtom atom1 = sortedAtomsByAtomCode.get(i);
				if (atom1.getProperty("position") != null)
					continue;
				conflictedAtoms.add(atom1);
				String atomCode1 = atom1.getProperty("code2");
				for (int j = i+1; j < sortedAtomsByAtomCode.size(); j++) {
					IAtom atom2 = sortedAtomsByAtomCode.get(j);
					String atomCode2 = atom2.getProperty("code2");
					if (atomCode1.equals(atomCode2))
						conflictedAtoms.add(atom2);
					else
						break;
				}
				int nbOfConflicts = conflictedAtoms.size();
				conflictSolver(pseudoMol, depth, position);
				position += nbOfConflicts;
			}
		}
		return finish(codeLayers, numberOfRepetitions);
	}
	
	/**
	 * Solve conflicts between atom sharing the same atom code
	 * @param pseudoMol
	 * @param depth
	 * @param position
	 * @return
	 */
	private boolean conflictSolver(IAtomContainer pseudoMol, int depth, int position) {
		int size = conflictedAtoms.size();
		if (size > 1) {
			if (depth != 0) {
				//use previous layer
				if (previousLayerSolver(position))
					return true;
				else
					removeSolvedConflictedAtom();
			}
			else {
				//initialize position
				int tempPos = position + size - 1;
				for (IAtom a : conflictedAtoms) {
					a.setProperty("position", tempPos);
				}
			}
			if (connectedBondsSolver())
				return true;
			else
				removeSolvedConflictedAtom();
			if (currentLayerSolver(position, pseudoMol))
				return true;
			else
				removeSolvedConflictedAtom();
			if (nextLayersSolver(pseudoMol))
				return true;
			//Symmetry, it does not matter. The order will not change the position. A position is attributed
			else {
				removeSolvedConflictedAtom();
				forceSolvedConflictedAtom();
			}
		}
		else {
			conflictedAtoms.get(0).setProperty("position", position);
			return true;
		}
		return false;
	}
	
	/**
	 * 
	 */
	private void forceSolvedConflictedAtom() {
		List<IAtom> temp = new ArrayList<IAtom>();
		Collections.sort(conflictedAtoms, new CompareByPosition());
		int currentPosition = -1;
		for (IAtom a : conflictedAtoms) {
			if (currentPosition == -1) {
				currentPosition = (int)a.getProperty("position");
				temp.add(a);
			}
			else {
				int nextAtomPos = (int)a.getProperty("position");
				if (currentPosition == nextAtomPos) {
					temp.add(a);
				}
				//attribute position for all atoms in temp 
				if (currentPosition != nextAtomPos || temp.size() == conflictedAtoms.size()) {
					int size = temp.size();
					for (int i = size-1; i > -1; i--) {
						IAtom tatom = temp.get(i);
						int newPosition = (int) tatom.getProperty("position") - (size - i - 1);
						tatom.setProperty("hasUniquePosition",true);
						tatom.setProperty("position",newPosition);
					}
					temp.clear();
				}
				else   {
					temp.add(a);
				}
				currentPosition = nextAtomPos;
			}
		}
		conflictedAtoms.clear();
	}
	
	/**
	 * Attribute a unique position
	 */
	private void removeSolvedConflictedAtom() {
		List<IAtom> newConflictedAtoms = new ArrayList<IAtom>();
		for (IAtom a : conflictedAtoms) {
			if ((boolean)a.getProperty("hasUniquePosition") == false) {
				newConflictedAtoms.add(a);
			}
		}
		conflictedAtoms = newConflictedAtoms;
	}
	
	/**
	 * @param position
	 * @return
	 */
	private boolean previousLayerSolver(int position) {
		boolean success = true;
		//sorted connection in previous layers
		for (int i = 0; i < conflictedAtoms.size(); i++) {
			int position2 = position + conflictedAtoms.size() - 1;
			IAtom atom1 = conflictedAtoms.get(i);
			atom1.setProperty("hasUniquePosition", true);
			List<IAtom> connectedAtomInPreviousLayerList1 = atom1.getProperty("connectedAtomInPreviousLayerList");
			List<Integer> priorities1 = new ArrayList<Integer>();
			for (IAtom a : connectedAtomInPreviousLayerList1) {
				priorities1.add(a.getProperty("position"));
			}
			Collections.sort(priorities1);
			for (int j = 0; j < conflictedAtoms.size(); j++) {
				if (i == j) {
					continue;
				}
				IAtom atom2 = conflictedAtoms.get(j);
				List<Integer> priorities2;
				if (atom2.getProperty("priorities") != null) {
					priorities2 = atom2.getProperty("priorities");
				}
				else {
					priorities2 = new ArrayList<Integer>();
					List<IAtom> connectedAtomInPreviousLayerList2 = atom2.getProperty("connectedAtomInPreviousLayerList");
					for (IAtom a : connectedAtomInPreviousLayerList2) {
						priorities2.add(a.getProperty("position"));
					}
					Collections.sort(priorities2);
					atom2.setProperty("priorities",priorities2);
				}
				
				//assign a priority based on priorities lists
				//connected to the same atoms in the previous layer
				if (priorities1.equals(priorities2)) {
					success = false;
					atom1.setProperty("hasUniquePosition", false);
				}
				else {
					if (priorities2.size() > priorities1.size()) {
						for (int a = 0; a < priorities2.size(); a++) {
							if (a < priorities1.size()){
								if (priorities1.get(a) < priorities2.get(a)) {
									position2--;
									break;
								}
							}
						}
					}
					else {
						for (int a = 0; a < priorities1.size(); a++) {
							if (a < priorities2.size()) {
								if (priorities1.get(a) < priorities2.get(a)) {
									position2--;
									break;
								}
							}
							//it means that atom 1 has one more connection with the previous layer and other connections are the same than atom 2
							else {
								position2--;
								break;
							}
						}
					}
				}
			}
			atom1.setProperty("position", position2);
		}
		return success;
	}

	/**
	 * @param position
	 * @param pseudoMol
	 * @return
	 */
	private boolean currentLayerSolver(int position, IAtomContainer pseudoMol) {
		boolean success = true;
		//sorted connection in previous layers
		for (int i = 0; i < conflictedAtoms.size(); i++) {
			int position2 = position + conflictedAtoms.size() - 1;
			IAtom atom1 = conflictedAtoms.get(i);
			atom1.setProperty("hasUniquePosition", true);
			int currentDepth = atom1.getProperty("depth");
			List<IAtom> connectedAtomPresentInCurrentLayer1 = 
					getConnectedAtomPresentInCurrentLayer(currentDepth, pseudoMol.getConnectedAtomsList(atom1));
			List<Integer> priorities1 = new ArrayList<Integer>();
			for (IAtom a : connectedAtomPresentInCurrentLayer1) {
				if (a.getProperty("position") != null)
					priorities1.add(a.getProperty("position"));
			}
			Collections.sort(priorities1);
			for (int j = 0; j < conflictedAtoms.size(); j++) {
				if (i == j) {
					continue;
				}
				IAtom atom2 = conflictedAtoms.get(j);
				List<Integer> priorities2;
				if (atom2.getProperty("priorities") != null) {
					priorities2 = atom2.getProperty("priorities");
				}
				else {
					priorities2 = new ArrayList<Integer>();
					List<IAtom> connectedAtomPresentInCurrentLayer2 = 
							getConnectedAtomPresentInCurrentLayer(currentDepth, pseudoMol.getConnectedAtomsList(atom2));
					for (IAtom a : connectedAtomPresentInCurrentLayer2) {
						if (a.getProperty("position") != null)
							priorities2.add(a.getProperty("position"));
					}
					Collections.sort(priorities2);
					atom2.setProperty("priorities",priorities2);
				}

				//assign a priority based on priorities lists
				//connected to the same atoms in the previous layer
				if (priorities1.equals(priorities2)) {
					success = false;
					atom1.setProperty("hasUniquePosition", false);
				}
				else {
					if (priorities2.size() > priorities1.size()) {
						for (int a = 0; a < priorities2.size(); a++) {
							if (a < priorities1.size()){
								if (priorities1.get(a) < priorities2.get(a)) {
									position2--;
									break;
								}
							}
						}
					}
					else {
						for (int a = 0; a < priorities1.size(); a++) {
							if (a < priorities2.size()) {
								if (priorities1.get(a) < priorities2.get(a)) {
									position2--;
									break;
								}
							}
							//it means that atom 1 has one more connection with the previous layer and other connections are the same than atom 2
							else {
								position2--;
								break;
							}
						}
					}
				}
			}
			atom1.setProperty("position", position2);
		}
		return success;
	}
	
	/**
	 * @return
	 */
	private boolean connectedBondsSolver() {
		boolean success = true;
		List<IAtom> sortedAtomsByABCode = new ArrayList<IAtom>(conflictedAtoms);
		//Normally in sortedAtomsByABCode the first elt has the lowest position property and the last the higher position property
		Collections.sort(sortedAtomsByABCode, new CompareByABCode());
		int size = sortedAtomsByABCode.size();
		//look by bonds codes
		for (int i = 0; i < size; i++) {
			IAtom atom1 = sortedAtomsByABCode.get(i);
			String abCode1 = atom1.getProperty("abCode");
			int currentPosition = (int) atom1.getProperty("position");
			int newPosition = (int) atom1.getProperty("position");
			atom1.setProperty("hasUniquePosition", true);
			for (int j = 0; j < size; j++) {
				if (i == j) 
					continue;
				IAtom atom2 = sortedAtomsByABCode.get(j);
				String abCode2 = atom2.getProperty("abCode");
				if (abCode1.equals(abCode2)) {
					success = false;
					atom1.setProperty("hasUniquePosition", false);
				}
				else {
					if (currentPosition == (int) atom2.getProperty("position"))
						newPosition--;
				}
			}
			atom1.setProperty("position", newPosition);
		}
		return success;
	}
	
	/**
	 * @param pseudoMol
	 * @return
	 */
	private boolean nextLayersSolver(IAtomContainer pseudoMol) {
		boolean success = true;
		int size = conflictedAtoms.size();
		//compare abCode of connected atoms in the next layers
		for (int i = 0; i < size; i++) {
			IAtom atom1 = conflictedAtoms.get(i);
			atom1.setProperty("hasUniquePosition", true);
			int currentDepth = atom1.getProperty("depth");
			List<IAtom> connectedAtoms1 = getConnectedAtomPresentInNextLayer(currentDepth, 
					pseudoMol.getConnectedAtomsList(atom1));
			Collections.sort(connectedAtoms1, new CompareByABCode());
			String abCode1 = "";
			for (IAtom a : connectedAtoms1) {
				abCode1 += a.getProperty("abCode");
			}
			int currentPosition = (int) atom1.getProperty("position");
			int newPosition = (int) atom1.getProperty("position");
			for (int j = 0; j < size; j++) {
				if (i == j) 
					continue;
				IAtom atom2 = conflictedAtoms.get(j);
				if (currentPosition == (int) atom2.getProperty("position")) {
					List<IAtom> connectedAtoms2 = getConnectedAtomPresentInNextLayer(atom2.getProperty("depth"), 
							pseudoMol.getConnectedAtomsList(atom2));
					Collections.sort(connectedAtoms2, new CompareByABCode());
					String abCode2 = "";
					for (IAtom a : connectedAtoms2) {
						abCode2 += a.getProperty("abCode");
					}
					if (abCode1.equals(abCode2)) {
						int res = nextLayersSolver(pseudoMol, connectedAtoms1, connectedAtoms2, currentDepth);
						//-9 is used when the next layers can't attribute a position (it's usually a symmetry)
						if (res == -9) {
							atom1.setProperty("hasUniquePosition", false);
							success = false;
						}
						else {
							newPosition = newPosition + res;
						}
					}
					else {
						if (currentPosition == (int) atom2.getProperty("position"))
							newPosition--;
					}
				}
				
			}
			atom1.setProperty("position", newPosition);
		}
		return success;
	}
	
	/**
	 * @param pseudoMol
	 * @param l1
	 * @param l2
	 * @param depth
	 * @return
	 */
	private int nextLayersSolver(IAtomContainer pseudoMol, List<IAtom> l1, List<IAtom> l2, int depth) {
		List<String> abCodes1 = new ArrayList<String>();
		List<String> abCodes2= new ArrayList<String>();
		List<IAtom> connectedAtoms1 = new ArrayList<IAtom>();
		List<IAtom> connectedAtoms2 = new ArrayList<IAtom>();
		for (IAtom a1 : l1) {
			connectedAtoms1 = getConnectedAtomPresentInNextLayer(depth, 
					pseudoMol.getConnectedAtomsList(a1));
			for (IAtom ca : connectedAtoms1) {
				abCodes1.add(ca.getProperty("abCode"));
			}
		}
		for (IAtom a2 : l2) {
			connectedAtoms2 = getConnectedAtomPresentInNextLayer(depth, 
					pseudoMol.getConnectedAtomsList(a2));
			for (IAtom ca : connectedAtoms2) {
				abCodes2.add(ca.getProperty("abCode"));
			}
		}
		Collections.sort(abCodes1, Collections.reverseOrder());
		Collections.sort(abCodes2, Collections.reverseOrder());
		if (abCodes1.equals(abCodes2) && !abCodes1.isEmpty()) {
			return nextLayersSolver(pseudoMol, connectedAtoms1, connectedAtoms2, depth+1); 
		}
		else {
			if (abCodes1.isEmpty() && !abCodes2.isEmpty()) {
				return 0;
			}
			else if (!abCodes1.isEmpty() && abCodes2.isEmpty()) {
				return -1;
			}
			else if (abCodes1.isEmpty() && abCodes2.isEmpty()) {
				return -9;
			}
			else {
				if (abCodes2.size() > abCodes1.size()) {
					for (int i = 0; i < abCodes1.size(); i++) {
						String abCode1 = abCodes1.get(i);
						String abCode2 = abCodes2.get(i);
						int compare = abCode1.compareTo(abCode2);
						if (compare == 0)
							continue;
						else if (compare > 0)
							return -1;
						else
							return 0;
					}
					return 0;
				}
				else {
					for (int i = 0; i < abCodes2.size(); i++) {
						String abCode1 = abCodes1.get(i);
						String abCode2 = abCodes2.get(i);
						int compare = abCode1.compareTo(abCode2);
						if (compare == 0)
							continue;
						else if (compare > 0)
							return -1;
						else
							return 0;
					}
					return -1;
				}
			}
		}
	}
	
	/**
	 * @param depth
	 * @param connectedAtom
	 * @return
	 */
	private List<IAtom> getConnectedAtomPresentInCurrentLayer(int depth, List<IAtom> connectedAtom) {
		List<IAtom> connectedAtomsInNextLayer = new ArrayList<IAtom>();
		for (IAtom a : connectedAtom) {
			if ((int)a.getProperty("depth") == depth) 
				connectedAtomsInNextLayer.add(a);
		}
		return connectedAtomsInNextLayer;
	}
	
	/**
	 * @param depth
	 * @param connectedAtom
	 * @return
	 */
	private List<IAtom> getConnectedAtomPresentInNextLayer(int depth, List<IAtom> connectedAtom) {
		List<IAtom> connectedAtomsInNextLayer = new ArrayList<IAtom>();
		for (IAtom a : connectedAtom) {
			if ((int)a.getProperty("depth") > depth) 
				connectedAtomsInNextLayer.add(a);
		}
		return connectedAtomsInNextLayer;
	}
	
	
	/**
	 * @param sphereList
	 * @param numberOfRepetitions
	 * @return
	 */
	private Map<String,String> finish(Map<Integer, Set<IAtom>> sphereList, Map<String, Integer> numberOfRepetitions) {
		Map<Integer,String> tempResult = new TreeMap<Integer,String>();
		int leavingIndicator = sphereList.size()*2;
		
		Map<IAtom, Integer> positionsOfStayingAtoms = new HashMap<IAtom, Integer>();
		Map<IAtom, Integer> positionsOfLeavingAtoms = new HashMap<IAtom, Integer>();
		
		for (Entry<Integer, Set<IAtom>> e : sphereList.entrySet()) {
			StringBuilder stayingCodes = new StringBuilder();
			StringBuilder leavingCodes = new StringBuilder();
			String[] stayingComplementCodes = new String[4];
			String[] leavingComplementCodes = new String[4];
			
			//initialize
			stayingComplementCodes[0] = "";
			stayingComplementCodes[1] = "";
			stayingComplementCodes[2] = "";
			stayingComplementCodes[3] = "";
			leavingComplementCodes[0] = "";
			leavingComplementCodes[1] = "";
			leavingComplementCodes[2] = "";
			leavingComplementCodes[3] = "";
			
			int depth = e.getKey();
			List<IAtom> atoms = new ArrayList<IAtom>(e.getValue());
			Collections.sort(atoms, new CompareByPosition());
			List<IBond> bondToAdd = new ArrayList<IBond>();
			
			int cptStaying = 0;
			int cptLeaving = 0;
			for (IAtom atom : atoms) {
				if (depth > 0) {
					if (atom.getProperty(additionalConstants.LEAVING_ATOM) != null) {
						positionsOfLeavingAtoms.put(atom, positionsOfLeavingAtoms.size());
					}
					else {
						positionsOfStayingAtoms.put(atom, positionsOfStayingAtoms.size());
					}
				}
				else {
					int pos = positionsOfStayingAtoms.size();
					positionsOfStayingAtoms.put(atom, pos);
					if (atom.getProperty("notInProduct") != null)
						atomInReactionCenterNotInProducts.add(pos+1);
				}
				List<IBond> bonds = mol.getConnectedBondsList(atom);
				for (IBond b : bonds) {
					IAtom other = b.getOther(atom);
					if (positionsOfStayingAtoms.containsKey(other) || positionsOfLeavingAtoms.containsKey(other)) {
						bondToAdd.add(b);
					}
				}
				StringBuilder code = new StringBuilder();
				//add atom
				String[] atomCode = atom.getProperty("code").toString().split("_")[0].split("/");

				String atomComplement = ""; 
				if (atomCode.length > 1) {
					atomComplement = atomCode[1];
				}

				code.append(atomCode[0]);
				code.append("(");
				
				//add bonds
				List<String> bondCodes = new ArrayList<String>();
				for (IBond b : bondToAdd) {
					String bondCode = b.getProperty("code").toString().split("_")[0] + "|";
					String pos = "";
					if (positionsOfStayingAtoms.get(b.getOther(atom)) != null) {
						pos = connectionTableAlphabet.get(positionsOfStayingAtoms.get(b.getOther(atom)));
					}
					else if (positionsOfLeavingAtoms.get(b.getOther(atom)) != null) {
						pos = String.format("%02X", positionsOfLeavingAtoms.get(b.getOther(atom)));
					}
					if (pos.length() > 0) {
						bondCode += pos;
						bondCodes.add(bondCode);
					}
					
				}
				Collections.sort(bondCodes, Collections.reverseOrder());
				int cptBond = 0;
				Map<Integer,String> bondComplements = new HashMap<Integer,String>();
				for (String bCode : bondCodes) {
					String[] codeAndConnections = bCode.split("\\|");
					String bondCode = codeAndConnections[0];
					String connections = "";
					if (codeAndConnections.length > 1) {
						connections = codeAndConnections[1];
					}
					String[] bondCodeWithComp = bondCode.split("/");
					bondCode = bondCodeWithComp[0];
					if (bondType == false) {
						bondCode = "";
					}
					bondCode += connections;
					if (bondCodeWithComp.length > 1) {
						bondComplements.put(cptBond, bondCodeWithComp[1]);
					}

					code.append(bondCode);
					cptBond++;
				}
				code.append(")");
				if (repetition == true) {
					code.append("["+numberOfRepetitions.get(atom.getID())+"]");
				}
				
				if (depth > 0) {
					if (atom.getProperty(additionalConstants.LEAVING_ATOM) != null) {
						leavingCodes.append(code.toString());
						
						//add complementary information (stereo and charge)
						String charges = leavingComplementCodes[0];
						String stereo = leavingComplementCodes[1];
						String isotopes = leavingComplementCodes[2];
						String radicals = leavingComplementCodes[3];
						if (atomComplement.length() > 0) {
							if (atomComplement.contains("~")) {
								int indexSymbol = atomComplement.indexOf("~");
								radicals += String.format("%02d", cptStaying) + atomComplement.substring(indexSymbol+1, indexSymbol+3);
								atomComplement = atomComplement.replace(atomComplement.substring(indexSymbol, indexSymbol+3),"");
							}
							if (atomComplement.contains("#")) {
								int indexSymbol = atomComplement.indexOf("#");
								isotopes += String.format("%02d", cptStaying) + atomComplement.substring(indexSymbol+1, indexSymbol+3);
								atomComplement = atomComplement.substring(0, indexSymbol);
							}
							if (charge == true && stereochemistry == true && atomComplement.length() > 0) {
								if (atomComplement.contains("@")) {
									if (atomComplement.indexOf("@") > 0) {
										charges += String.format("%02d", cptLeaving) + atomComplement.substring(0, 2);
										stereo += String.format("%02d", cptLeaving) + atomComplement.substring(3);
									}
									else {
										stereo += String.format("%02d", cptLeaving) + atomComplement.substring(1);
									}
								}
								else {
									charges += String.format("%02d", cptLeaving) + atomComplement.substring(0);

								}
							}
							else if (charge == true && stereochemistry == false && atomComplement.length() > 0) {
								charges += String.format("%02d", cptLeaving) + atomComplement;
							}
							else if (charge == false && stereochemistry == true && atomComplement.length() > 0) {
								stereo += String.format("%02d", cptLeaving) + atomComplement;
							}
						}
						//increment by 1 because 1 atom is added
						cptLeaving++;
						for (Entry<Integer,String> e2 : bondComplements.entrySet()) {
							int index = e2.getKey() + cptLeaving;
							stereo += String.format("%02d", index) + e2.getValue();
						}
						leavingComplementCodes[0] = charges;
						leavingComplementCodes[1] = stereo;
						leavingComplementCodes[2] = isotopes;
						leavingComplementCodes[3] = radicals;
						//increment with the number of bonds (corresponds to the position of the next atom)
						cptLeaving = cptLeaving + cptBond;
					}
					else {
						stayingCodes.append(code.toString());
						
						//add complementary information (stereo and charge)
						String charges = stayingComplementCodes[0];
						String stereo = stayingComplementCodes[1];
						String isotopes = stayingComplementCodes[2];
						String radicals = stayingComplementCodes[3];
						if (atomComplement.length() > 0) {
							if (atomComplement.contains("~")) {
								int indexSymbol = atomComplement.indexOf("~");
								radicals += String.format("%02d", cptStaying) + atomComplement.substring(indexSymbol+1, indexSymbol+3);
								atomComplement = atomComplement.replace(atomComplement.substring(indexSymbol, indexSymbol+3),"");
							}
							if (atomComplement.contains("#")) {
								int indexSymbol = atomComplement.indexOf("#");
								isotopes += String.format("%02d", cptStaying) + atomComplement.substring(indexSymbol+1, indexSymbol+3);
								atomComplement = atomComplement.substring(0, indexSymbol);
							}
							if (charge == true && stereochemistry == true && atomComplement.length() > 0) {
								if (atomComplement.contains("@")) {
									if (atomComplement.indexOf("@") > 0) {
										charges += String.format("%02d", cptStaying) + atomComplement.substring(0, 2);
										stereo += String.format("%02d", cptStaying) + atomComplement.substring(3);
									}
									else {
										stereo += String.format("%02d", cptStaying) + atomComplement.substring(1);
									}
								}
								else {
									charges +=String.format("%02d",  cptStaying) + atomComplement.substring(0);

								}
							}
							else if (charge == true && stereochemistry == false && atomComplement.length() > 0) {
								charges += String.format("%02d", cptStaying) + atomComplement;
							}
							else if  (charge == false && stereochemistry == true && atomComplement.length() > 0) {
								stereo += String.format("%02d", cptStaying) + atomComplement;
							}
						}
						//increment by 1 because 1 atom is added
						cptStaying++;
						for (Entry<Integer,String> e2 : bondComplements.entrySet()) {
							int index = e2.getKey() + cptStaying;
							stereo += String.format("%02d", index) + e2.getValue();
						}
						stayingComplementCodes[0] = charges;
						stayingComplementCodes[1] = stereo;
						stayingComplementCodes[2] = isotopes;
						stayingComplementCodes[3] = radicals;
						//increment with the number of bonds (corresponds to the position of the next atom)
						cptStaying = cptStaying + cptBond;
					}
				}
				else {
					stayingCodes.append(code.toString());
					
					//add complementary information (stereo and charge)
					String charges = stayingComplementCodes[0];
					String stereo = stayingComplementCodes[1];
					String isotopes = stayingComplementCodes[2];
					String radicals = stayingComplementCodes[3];
					if (atomComplement.length() > 0) {
						if (atomComplement.contains("~")) {
							int indexSymbol = atomComplement.indexOf("~");
							radicals += String.format("%02d", cptStaying) + atomComplement.substring(indexSymbol+1, indexSymbol+3);
							atomComplement = atomComplement.replace(atomComplement.substring(indexSymbol, indexSymbol+3),"");
						}
						if (atomComplement.contains("#")) {
							int indexSymbol = atomComplement.indexOf("#");
							isotopes += String.format("%02d", cptStaying) + atomComplement.substring(indexSymbol+1, indexSymbol+3);
							atomComplement = atomComplement.substring(0, indexSymbol);
						}
						if (charge == true && stereochemistry == true && atomComplement.length() > 0) {
							if (atomComplement.contains("@")) {
								if (atomComplement.indexOf("@") > 0) {
									charges +=String.format("%02d",  cptStaying) + atomComplement.substring(0, 2);
									stereo += String.format("%02d", cptStaying) + atomComplement.substring(3);
								}
								else {
									stereo += String.format("%02d", cptStaying) + atomComplement.substring(3);
								}
							}
							else {
								charges += String.format("%02d", cptStaying) + atomComplement.substring(0);

							}
						}
						else if (charge == true && stereochemistry == false && atomComplement.length() > 0) {
							charges += String.format("%02d", cptStaying) + atomComplement;
						}
						else if  (charge == false && stereochemistry == true && atomComplement.length() > 0) {
							stereo += String.format("%02d", cptStaying) + atomComplement.replace("@","");
						}
						
					}
					//increment by 1 because 1 atom is added
					cptStaying++;
					for (Entry<Integer,String> e2 : bondComplements.entrySet()) {
						int index = e2.getKey() + cptStaying;
						stereo += String.format("%02d", index) + e2.getValue();
					}
					stayingComplementCodes[0] = charges;
					stayingComplementCodes[1] = stereo;
					stayingComplementCodes[2] = isotopes;
					stayingComplementCodes[3] = radicals;
					//increment with the number of bonds (corresponds to the position of the next atom)
					cptStaying = cptStaying + cptBond;
				}
			}
			if (stayingCodes.length() > 0) {
				if (stayingComplementCodes[0].length() > 0) {
					stayingCodes.append("/c"+stayingComplementCodes[0]);
				}
				if (stayingComplementCodes[1].length() > 0) {
					stayingCodes.append("/s"+stayingComplementCodes[1]);
				}
				if (stayingComplementCodes[2].length() > 0) {
					stayingCodes.append("/i"+stayingComplementCodes[2]);
				}
				if (stayingComplementCodes[3].length() > 0) {
					stayingCodes.append("/r"+stayingComplementCodes[3]);
				}
				tempResult.put(depth, stayingCodes.toString());
			}
			if (leavingCodes.length() > 0) {
				if (leavingComplementCodes[0].length() > 0) {
					leavingCodes.append("/c"+leavingComplementCodes[0]);
				}
				if (leavingComplementCodes[1].length() > 0) {
					leavingCodes.append("/s"+leavingComplementCodes[1]);
				}
				if (leavingComplementCodes[2].length() > 0) {
					leavingCodes.append("/i"+leavingComplementCodes[2]);
				}
				if (leavingComplementCodes[3].length() > 0) {
					leavingCodes.append("/r"+leavingComplementCodes[3]);
				}
				tempResult.put(leavingIndicator+depth, leavingCodes.toString());

			}
		}
		
		//sort result by depth and group
		Map<String,String> result = new LinkedHashMap<String,String>();
		for (Entry<Integer,String> e : tempResult.entrySet()) {
			int depth2 = e.getKey();
			if (depth2 < leavingIndicator) {
				result.put(String.valueOf(depth2), e.getValue());
			}
			else {
				int newDepth = depth2 - leavingIndicator;
				result.put(Character.toString((char) (64+newDepth)), e.getValue());
			}
		}
		return result;
	}
	

	/**
	 * @param acSet
	 * @return
	 */
	private HashMap<String, IAtom> indexAtomIAtomContainerSet (IAtomContainerSet acSet) {
		HashMap<String, IAtom> result = new HashMap<String, IAtom>();

		for (int i = 0; i < acSet.getAtomContainerCount(); i ++) {
			IAtomContainer ac = acSet.getAtomContainer(i);
			if ((boolean)ac.getProperty(additionalConstants.SPECTATOR) == true)
				continue;
			for (int j = 0; j < ac.getAtomCount(); j++) {
				result.put(ac.getAtom(j).getID(), ac.getAtom(j));
			}
		}
		return result;
	}
	
	/**
	 * @param acSet
	 * @return
	 */
	private HashMap<String, IBond> indexAtomIBondContainerSet(IAtomContainerSet acSet) {
		HashMap<String, IBond> result = new HashMap<String, IBond>();

		for (int i = 0; i < acSet.getAtomContainerCount(); i ++) {
			IAtomContainer ac = acSet.getAtomContainer(i);
			for (int j = 0; j < ac.getBondCount(); j++) {
				result.put(ac.getBond(j).getID(), ac.getBond(j));
			}
		}
		return result;
	}

	/**
	 * @param atomInProduct
	 * @param bonds
	 * @param indexAtomReactants
	 * @return
	 */
	private String encodeAtom(IAtom atomInProduct, List<IBond> bonds, Map<String, IAtom> indexAtomReactants) {
		StringBuilder code = new StringBuilder();
		String bScore = Integer.toHexString(getRCScore(bonds));
		code.append(bScore);
		IAtom atomInReactant = indexAtomReactants.get(atomInProduct.getProperty(ATOM_ATOM_MAPPING)+"");
		
		if (hybridization == true) {
			String rHybridization = "";
			String pHybridization = "";
			if (atomInReactant != null) {
				rHybridization = hybridizationEncoding(atomInReactant.getHybridization());
				pHybridization = hybridizationEncoding(atomInProduct.getHybridization());
			}
			else {
				rHybridization = "0";
				pHybridization = hybridizationEncoding(atomInProduct.getHybridization());
			}
			code.append(rHybridization);
			code.append(pHybridization);
		}
		
		String encSymbol = encodeSymbol(atomInProduct);
		code.append(encSymbol);
		
		code.append("/");
		if (charge == true) {
			String chargeInReactant = "";
			String chargeInProduct = "";
			if (atomInReactant != null) {
				chargeInReactant = encodeChargeRadicalAndIsotope(atomInReactant.getFormalCharge());
				chargeInProduct = encodeChargeRadicalAndIsotope(atomInProduct.getFormalCharge());
			}
			else {
				chargeInReactant = "0";
				chargeInProduct = encodeChargeRadicalAndIsotope(atomInProduct.getFormalCharge());
			}
			if (!chargeInReactant.equals("0") || !chargeInProduct.equals("0")) {
				code.append(chargeInReactant);
				code.append(chargeInProduct);
			}
		}
		if (stereochemistry = true) {
			String stereoInReactant = encodeAtomStereo(atomInReactant);
			String stereoInProduct = encodeAtomStereo(atomInProduct);
			if (!stereoInReactant.equals("0") || !stereoInProduct.equals("0")) {
				code.append("@");
				code.append(stereoInReactant);
				code.append(stereoInProduct);
			}
			
		}
		
		String radicalInReactant = "0";
		String radicalInProduct = "0";
		if (atomInProduct.getProperty(additionalConstants.RADICAL) != null) {
			int radical = atomInProduct.getProperty(additionalConstants.RADICAL);
			radicalInProduct = encodeChargeRadicalAndIsotope(radical);
		}
		if (atomInReactant.getProperty(additionalConstants.RADICAL) != null) {
			int radical = atomInReactant.getProperty(additionalConstants.RADICAL);
			radicalInReactant = encodeChargeRadicalAndIsotope(radical);
		}
		if (!radicalInReactant.equals("0") || !radicalInProduct.equals("0")) {
			code.append("~");
			code.append(radicalInReactant);
			code.append(radicalInProduct);
		}
		
		String isotopeInReactant = "0";
		String isotopeInProduct = "0";
		if (atomInProduct.getMassNumber() != null && !atomInProduct.getSymbol().equals("H")) {
			int frequentMass = ElementCalculation.calculateMass(atomInProduct.getSymbol());
			if (frequentMass != atomInProduct.getMassNumber()) {
				int delta = atomInProduct.getMassNumber() - frequentMass;
				isotopeInProduct = encodeChargeRadicalAndIsotope(delta);
			}
		}
		if (atomInReactant.getMassNumber() != null && !atomInReactant.getSymbol().equals("H")) {
			int frequentMass = ElementCalculation.calculateMass(atomInReactant.getSymbol());
			if (frequentMass != atomInReactant.getMassNumber()) {
				int delta = atomInReactant.getMassNumber() - frequentMass;
				isotopeInReactant = encodeChargeRadicalAndIsotope(delta);
			}
		}
		if (!isotopeInReactant.equals("0") || !isotopeInReactant.equals("0")) {
			code.append("#");
			code.append(isotopeInReactant);
			code.append(isotopeInProduct);
		}
		//will be use just to sort atom
		code.append("_");
		code.append(encodeBoolean(atomInReactant.isInRing()));
		code.append(encodeBoolean(atomInReactant.isAromatic()));
		code.append(encodeBoolean(atomInProduct.isInRing()));
		code.append(encodeBoolean(atomInProduct.isAromatic()));
		return code.toString();
	}
	
	/**
	 * @param inProduct
	 * @param indexBondReactants
	 * @return
	 */
	private String encodeBond(IBond inProduct, Map<String, IBond> indexBondReactants) {
		StringBuilder code = new StringBuilder();
		IBond inReactant = indexBondReactants.get(inProduct.getID());
		
		if (inProduct.getProperty(additionalConstants.BOND_CHANGE_INFORMATION) != null) {
			if ((int) inProduct.getProperty(additionalConstants.BOND_CHANGE_INFORMATION) == additionalConstants.BOND_MADE) {
				inReactant = new Bond();
				inReactant.setIsAromatic(false);
				inReactant.setProperty("ignore", true);
			}
		}
		
		if (inReactant.getProperty(additionalConstants.BOND_CHANGE_INFORMATION) != null) {
			if ((int) inReactant.getProperty(additionalConstants.BOND_CHANGE_INFORMATION) == additionalConstants.BOND_CLEAVED) {
				inProduct = new Bond();
				inProduct.setIsAromatic(false);
				inProduct.setProperty("ignore", true);
			}
		}
		
		if (inReactant.isAromatic()) {
			code.append("9");
		}
		else {
			//bond formed (not in reactant)
			if (inReactant.getProperty("ignore") != null) {
				code.append("0");
			}
			else {
				code.append(inReactant.getOrder().numeric());
			}
		}

		if (inProduct.isAromatic()) {
			code.append("9");
		}
		else {
			//bond cleaved (not in product)
			if (inProduct.getProperty("ignore") != null) {
				code.append("0");
			}
			else {
				code.append(inProduct.getOrder().numeric());
			}
		}
		
		code.append("/");
		if (stereochemistry == true) {
			String stereoInReactant;
			String stereoInProduct;
			if (inReactant.getProperty("ignore") != null) {
				stereoInReactant = "0";
			}
			else {
				stereoInReactant = encodeBondStereo(inReactant);
			}
			
			if (inProduct.getProperty("ignore") != null) {
				stereoInProduct = "0";
			}
			else {
				stereoInProduct = encodeBondStereo(inProduct);
			}

			if (!stereoInReactant.equals("0") || !stereoInProduct.equals("0")) {
				code.append(stereoInReactant);
				code.append(stereoInProduct);
			}
		}

		//will be use just to sort bonds
		code.append("_");
		code.append(encodeBoolean(inReactant.isInRing()));
		code.append(encodeBoolean(inReactant.isAromatic()));
		code.append(encodeBoolean(inProduct.isInRing()));
		code.append(encodeBoolean(inProduct.isAromatic()));
		
		return code.toString();
	}
	
	
	/**
	 * @param bond
	 * @return
	 */
	private String encodeBondStereo(IBond bond) {
		if (bond.getStereo().equals(IBond.Stereo.UP_OR_DOWN_INVERTED)) 
			return "9";
		else if (bond.getStereo().equals(IBond.Stereo.UP_OR_DOWN)) 
			return "8";
		else if (bond.getStereo().equals(IBond.Stereo.E_OR_Z)) 
			return "7";
		else if (bond.getStereo().equals(IBond.Stereo.UP)) 
			return "6";
		else if (bond.getStereo().equals(IBond.Stereo.DOWN)) 
			return "4";
		else if (bond.getStereo().equals(IBond.Stereo.UP_INVERTED)) 
			return "5";
		else if (bond.getStereo().equals(IBond.Stereo.DOWN_INVERTED)) 
			return "3";
		else if (bond.getStereo().equals(IBond.Stereo.Z)) 
			return "2";
		else if (bond.getStereo().equals(IBond.Stereo.E)) 
			return "1";
		else
			return "0";
	}

	/**
	 * @param b
	 * @return
	 */
	private int encodeBoolean(boolean b) {
		return b ? 1 : 0;
	}


	/**
	 * @param chargeOrMass
	 * @return
	 */
	private String encodeChargeRadicalAndIsotope(int chargeOrMass) {
		String chargeEncoded = "0";
		if (chargeOrMass == -17) {
			chargeEncoded = "1";
		}
		else if (chargeOrMass == -16) {
			chargeEncoded = "2";
		}
		else if (chargeOrMass == -15) {
			chargeEncoded = "3";
		}
		else if (chargeOrMass == -14) {
			chargeEncoded = "4";
		}
		else if (chargeOrMass == -13) {
			chargeEncoded = "5";
		}
		else if (chargeOrMass == -12) {
			chargeEncoded = "6";
		}
		else if (chargeOrMass == -11) {
			chargeEncoded = "7";
		}
		else if (chargeOrMass == -10) {
			chargeEncoded = "8";
		}
		else if (chargeOrMass == -9) {
			chargeEncoded = "9";
		}
		else if (chargeOrMass == -8) {
			chargeEncoded = "A";
		}
		else if (chargeOrMass == -7) {
			chargeEncoded = "B";
		}
		else if (chargeOrMass == -6) {
			chargeEncoded = "C";
		}
		else if (chargeOrMass == -5) {
			chargeEncoded = "D";
		}
		else if (chargeOrMass == -4) {
			chargeEncoded = "E";
		}
		else if (chargeOrMass == -3) {
			chargeEncoded = "F";
		}
		else if (chargeOrMass == -2) {
			chargeEncoded = "G";
		}
		else if (chargeOrMass == -1) {
			chargeEncoded = "H";
		}
		else if (chargeOrMass == 1) {
			chargeEncoded = "I";
		}
		else if (chargeOrMass == 2) {
			chargeEncoded = "J";
		}
		else if (chargeOrMass == 3) {
			chargeEncoded = "K";
		}
		else if (chargeOrMass == 4) {
			chargeEncoded = "L";
		}
		else if (chargeOrMass == 5) {
			chargeEncoded = "M";
		}
		else if (chargeOrMass == 6) {
			chargeEncoded = "N";
		}
		else if (chargeOrMass == 7) {
			chargeEncoded = "O";
		}
		else if (chargeOrMass == 8) {
			chargeEncoded = "P";
		}
		else if (chargeOrMass == 9) {
			chargeEncoded = "Q";
		}
		else if (chargeOrMass == 10) {
			chargeEncoded = "R";
		}
		else if (chargeOrMass == 11) {
			chargeEncoded = "S";
		}
		else if (chargeOrMass == 12) {
			chargeEncoded = "T";
		}
		else if (chargeOrMass == 13) {
			chargeEncoded = "U";
		}
		else if (chargeOrMass == 14) {
			chargeEncoded = "V";
		}
		else if (chargeOrMass == 15) {
			chargeEncoded = "W";
		}
		else if (chargeOrMass == 16) {
			chargeEncoded = "X";
		}
		else if (chargeOrMass == 17) {
			chargeEncoded = "Y";
		}

		return chargeEncoded;
	}
	

	/**
	 * @param atom
	 * @return
	 */
	private String encodeAtomStereo(IAtom atom) {
		if (atom.getProperty(additionalConstants.STEREO_TYPE) != null) {
			String type = atom.getProperty(additionalConstants.STEREO_TYPE);
			String shorthand = atom.getProperty(additionalConstants.STEREO_SHORTHAND);
			if (type.equals("Tetrahedral") && shorthand.equals("ANTI_CLOCKWISE")) {
				return "1";
			}
			else if (type.equals("Tetrahedral") && shorthand.equals("CLOCKWISE")) {
				return "2";
			}
			else if (type.equals("ExtendedTetrahedral") && shorthand.equals("ANTI_CLOCKWISE")) {
				return "3";
			}
			else if (type.equals("ExtendedTetrahedral") && shorthand.equals("CLOCKWISE")) {
				return "4";
			}
			if (type.equals("DoubleBond") && shorthand.equals("ANTI_CLOCKWISE")) {
				return "5";
			}
			if (type.equals("DoubleBond") && shorthand.equals("CLOCKWISE")) {
				return "6";
			}
			if (type.equals("ExtendedCisTrans") && shorthand.equals("ANTI_CLOCKWISE")) {
				return "7";
			}
			if (type.equals("ExtendedCisTrans") && shorthand.equals("CLOCKWISE")) {
				return "8";
			}
			else if (type.equals("SquarePlanar1")) {
				return "9";
			}
			else if (type.equals("SquarePlanar2")) {
				return "A";
			}
			else if (type.equals("SquarePlanar3")) {
				return "B";
			}
			else if (type.equals("TrigonalBipyramidal") && shorthand.equals("ANTI_CLOCKWISE")) {
				return "C";
			}
			else if (type.equals("TrigonalBipyramidal") && shorthand.equals("CLOCKWISE")) {
				return "D";
			}
			else if (type.equals("Octahedral") && shorthand.equals("ANTI_CLOCKWISE")) {
				return "E";
			}
			else if (type.equals("Octahedral") && shorthand.equals("CLOCKWISE")) {
				return "F";
			}
			else if (type.equals("Atropisomeric") && shorthand.equals("ANTI_CLOCKWISE")) {
				return "G";
			}
			else if (type.equals("Atropisomeric") && shorthand.equals("CLOCKWISE")) {
				return "H";
			}
		}
		return "0";
	}

	/**
	 * Generate individual atom and bond code and also abCode (atom code), which combines
	 * the atomCode and the reverse sorted code of all connected bonds (will be used for ranking)
	 * @param sphere
	 * @param indexAtomReactants
	 * @param indexBondReactants
	 * @param atomContainer
	 * @param sphereList
	 * @param visited
	 * @param depth
	 * @param cpt
	 */
	private void generateAllAtomAndBondCodes(Set<IAtom> sphere, Map<String, IAtom> indexAtomReactants, 
			Map<String, IBond> indexBondReactants, IAtomContainer atomContainer, Map<Integer, Set<IAtom>> sphereList, 
			List<IAtom> visited, int depth, int cpt) {

		IAtom nextAtom;
		Set<IAtom> newSphere = new HashSet<IAtom>();

		visited.addAll(sphere);
		for (IAtom atom : sphere) {
			String abCode = "";
			if (atom.getProperty("depth") == null) {
				atom.setProperty("depth", depth);
			}
			
			//generate atom code
			List<IBond> bonds = atomContainer.getConnectedBondsList(atom);
			if (atom.getProperty("code") == null) {
				String code = encodeAtom(atom, bonds, indexAtomReactants);
				atom.setProperty("code", code+"_"+cpt);
				atom.setProperty("code2", code);
				atom.setProperty("connectedBondsList", bonds);
				abCode += code;
			}
			else {
				abCode = atom.getProperty("code2");
			}
			
			//generate bond code and other part of abCode
			List<String> connectedBondCodes = new ArrayList<String>();
			List<IAtom> connectedAtomInPreviousLayerList = new ArrayList<IAtom>();
			for (IBond bond : bonds) {
				String bondCode;
				if (bond.getProperty("code") == null) {
					bondCode = encodeBond(bond, indexBondReactants);
					bond.setProperty("code", bondCode+"_"+depth);
					bond.setProperty("code2", bondCode);
				}
				else {
					bondCode = bond.getProperty("code2");
				}
				connectedBondCodes.add(bondCode);
				
				nextAtom = bond.getOther(atom);
				if (!visited.contains(nextAtom)) {
					if (!sphere.contains(nextAtom)) {
						newSphere.add(nextAtom);
					}
				}
				if (nextAtom.getProperty("depth") != null) {
					int d = (int) nextAtom.getProperty("depth");
					if (d == depth - 1) {
						connectedAtomInPreviousLayerList.add(nextAtom);
					}
				}
			}
			atom.setProperty("connectedAtomInPreviousLayerList", connectedAtomInPreviousLayerList);
			
			//complete abcode by reverse sorting bond code and merge all codes
			Collections.sort(connectedBondCodes, Collections.reverseOrder());
			for (String bCode : connectedBondCodes) {
				abCode += bCode;
			}
			atom.setProperty("abCode", abCode);
		}
		if (newSphere.size() > 0) {
			//dead code now because the depth start at 1 as atoms in RC are in the initial list
			/*
			if (depth == 0) {
				sphereList.put(0, newSphere);
			}
			*/
			sphereList.put(depth, newSphere);
			generateAllAtomAndBondCodes(newSphere, indexAtomReactants, indexBondReactants, atomContainer, 
					sphereList, visited, depth+1, cpt-1);
		}
	}

	/**
	 * @param hybridization
	 * @return
	 */
	private String hybridizationEncoding(IAtomType.Hybridization hybridization) {
		try {
			if (hybridization.equals(IAtomType.Hybridization.PLANAR3)) return Integer.toHexString(1);
			if (hybridization.equals(IAtomType.Hybridization.S)) return Integer.toHexString(2);
			if (hybridization.equals(IAtomType.Hybridization.SP1)) return Integer.toHexString(3);
			if (hybridization.equals(IAtomType.Hybridization.SP2)) return Integer.toHexString(4);
			if (hybridization.equals(IAtomType.Hybridization.SP3)) return Integer.toHexString(5);
			if (hybridization.equals(IAtomType.Hybridization.SP3D1)) return Integer.toHexString(6);
			if (hybridization.equals(IAtomType.Hybridization.SP3D2)) return Integer.toHexString(7);
			if (hybridization.equals(IAtomType.Hybridization.SP3D3)) return Integer.toHexString(8);
			if (hybridization.equals(IAtomType.Hybridization.SP3D4)) return Integer.toHexString(9);
			if (hybridization.equals(IAtomType.Hybridization.SP3D5)) return Integer.toHexString(10);
		}
		catch(NullPointerException e) {
			return Integer.toHexString(11);
		}
		return Integer.toHexString(11);
	}

	/**
	 * @param atom
	 * @return
	 */
	private String encodeSymbol(IAtom atom) {
		String code = null;
		String symbol = atom.getSymbol();
		if ("R".equals(symbol)) {
            code = "FF";
        } 
		else if ("*".equals(symbol)) {
            code = "FE";
        } 
		else if ("H".equals(symbol)) {
			int mass = 1;
			if (atom.getMassNumber() != null) {
				mass = atom.getMassNumber();
			}
			if (mass == 1) {
				code = String.format("%02X", PeriodicTable.getAtomicNumber(symbol));
			}
			//Deuterium
			else if (mass == 2) {
				code = "FD";
			}
			//Triterium
			else if (mass == 3) {
				code = "FC";
			}
        } 
        else {
        	code = String.format("%02X", PeriodicTable.getAtomicNumber(symbol));
        }
		return code;
	}
	
	/**
	 * 	Reacting center status
	 * 0 = unmarked, 1 = a center, -1 = not a center,
	 * Additional: 2 = no change,
	 * 4 = bond order changes
	 * 6 = bond broken,
	 * 8 = bond made,
	 * 10 = 4+8 (both broken and order changes);
	 * 12 = 6+8 (both made and order changes);
	 * 5 = (4 + 1), 7 = (6 + 1), 9 = (8 + 1), 11 = (10 + 1) and 13 = (12 + 1)
	 * NB: the score related to the higher bond change information is kept
	 * @param bond
	 * @return
	 */
	private int getRCScore(List<IBond> bonds) {
		int centerStatusValue = -1;
		for (IBond bond : bonds) {
			int score = 0;
			if (bond.getProperty(BOND_CHANGE_INFORMATION) != null) {
				score += 1;
				if ((int) bond.getProperty(BOND_CHANGE_INFORMATION) == additionalConstants.BOND_CLEAVED) 
					score += 6;
				if ((int) bond.getProperty(BOND_CHANGE_INFORMATION) == additionalConstants.BOND_MADE) 
					score += 8;
				if ((int) bond.getProperty(BOND_CHANGE_INFORMATION) == additionalConstants.BOND_ORDER) 
					score += 4;
			}
			else {
				score += 0;
			}
			if (score > 0 && score > centerStatusValue) {
				centerStatusValue = score;
			}
		}
		
		return (centerStatusValue != -1) ? centerStatusValue : 0;
	}
	

	/**
	 * @param map
	 * @return
	 */
	public String reactionCodeMapToStringOneLine(Map<String,String> map) {
		StringBuilder reactionCode = new StringBuilder();
		for (Entry<String,String> e : map.entrySet()) {
			reactionCode.append(e.getKey());
			reactionCode.append(":");
			reactionCode.append(e.getValue());
			reactionCode.append("|");
		}
		return reactionCode.toString();
	}
	
	/**
	 * @param map
	 * @return
	 */
	public String reactionCodeMapToStringMultiLines(Map<String,String> map) {
		StringBuilder reactionCode = new StringBuilder();
		for (Entry<String,String> e : map.entrySet()) {
			reactionCode.append(e.getKey());
			reactionCode.append(":");
			reactionCode.append(e.getValue());
			reactionCode.append("|\n");
		}
		return reactionCode.toString();
	}
	
	/**
	 * @param map
	 * @return
	 */
	/*private <K, V extends Comparable<? super V>> Map<K, V> sortByValue(Map<K, V> map) {
        List<Entry<K, V>> list = new ArrayList<>(map.entrySet());
        list.sort(Entry.comparingByValue());

        Map<K, V> result = new LinkedHashMap<>();
        for (Entry<K, V> entry : list) {
            result.put(entry.getKey(), entry.getValue());
        }

        return result;
    }
	*/
	
	
	/**
	 * @return
	 */
	private HashMap<Integer,String> makeConnectionTableAlphabet(){
		HashMap<Integer,String> result = new HashMap<Integer,String>();
		int cpt = 0;
		for (char alphabet1 = 'G'; alphabet1 <= 'Z'; alphabet1++) {
			for (char alphabet2 = 'G'; alphabet2 <= 'Z'; alphabet2++) {
				result.put(cpt, ""+alphabet1+alphabet2);
				cpt++;
			}
		}
		return result;
	}
	
	/**
	 * @param container
	 */
	private void reinitializedVisitedConstants(IAtomContainer container) {
		for (int i = 0; i < container.getAtomCount(); i++) {
			container.getAtom(i).setFlag(CDKConstants.VISITED, false);
		}
		for (int i = 0; i < container.getBondCount(); i++) {
			container.getBond(i).setFlag(CDKConstants.VISITED, false);
		}
	}
	
	/**
	 * @return
	 */
	public String getReactionCode() {
		return reactionCode;
	}
	
	/**
	 * @param hybridization
	 */
	public void setHybridization(boolean hybridization) {
		this.hybridization = hybridization;
	}
	
	/**
	 * @param stereochemistry
	 */
	public void setStereochemistry(boolean stereochemistry) {
		this.stereochemistry = stereochemistry;
	}
	
	/**
	 * @param charge
	 */
	public void setCharge(boolean charge) {
		this.charge = charge;
	}
	
	/**
	 * @param bondType
	 */
	public void setBondType(boolean bondType) {
		this.bondType = bondType;
	}
	
	/**
	 * @param repetition
	 */
	public void setRepetition(boolean repetition) {
		this.repetition = repetition;
	}
	
	/**
	 * @return
	 */
	public Set<Integer> getAtomInReactionCenterNotInProducts() {
		return atomInReactionCenterNotInProducts;
	}
	
}


//ascending comparison res = (1,2,3)
class CompareByPosition implements Comparator<IAtom> { 
	public int compare(IAtom a1, IAtom a2) { 
		if ((int) a1.getProperty("position") > (int) a2.getProperty("position"))
			return 1;
		else 
			return -1;
	} 
}

//descending comparison res = (3,2,1)
class CompareByAtomCode implements Comparator<IAtom> { 
	public int compare(IAtom a1, IAtom a2) { 
		if (((String) a1.getProperty("code")).compareTo(a2.getProperty("code")) < 0)
			return 1;
		else 
			return -1;
	} 
}

//descending comparison res = (3,2,1)
class CompareByABCode implements Comparator<IAtom> { 
	public int compare(IAtom a1, IAtom a2) { 
		if (((String) a1.getProperty("abCode")).compareTo(a2.getProperty("abCode")) < 0)
			return 1;
		else 
			return -1;
	} 
}
