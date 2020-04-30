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

import org.openscience.cdk.Bond;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
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

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}
	public EncodeReactionCode() {
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
	public Map<String,String> makeReactionCode(Set<IAtom> atomInRC, IAtomContainerSet reactants,
			IAtomContainer pseudoMolecule, Map<String, Integer> numberOfRepetitions) {

		mol = pseudoMolecule;
		reactionCode = null;
		
		//make table of correspondence to connect atom for reaction code generation
		connectionTableAlphabet = makeConnectionTableAlphabet();
		
		//map containing codes of each layer
		Map<String,String> reactionCodeMap = new HashMap<String,String>();
		
		//generate index
		Map<String, IAtom> indexAtomReactants = indexAtomIAtomContainerSet(reactants);
		Map<String, IBond> indexBondReactants = indexAtomIBondContainerSet(reactants);
		reinitializedConstants(pseudoMolecule);

		Map<Integer, Set<IAtom>> sphereList = new HashMap<Integer, Set<IAtom>>();
		sphereList.put(0, atomInRC);
		generateAllAtomAndBondCodes(atomInRC, indexAtomReactants, indexBondReactants, pseudoMolecule, sphereList, 
				new ArrayList<IAtom>(), 1, pseudoMolecule.getAtomCount());
		
		List<IAtom> sortedAtomInRC = new ArrayList<IAtom>(atomInRC);

		Collections.sort(sortedAtomInRC, new CompareByAtomCode());
		List<IAtom> roots = getRootAtoms(sortedAtomInRC);
		for (IAtom root : roots) {
			//reset Properties define to rank atoms
			resetRankingProperties();
			List<IAtom> sphere = new ArrayList<IAtom>();
//			List<IAtom> visited = new ArrayList<IAtom>();
//			for (Entry<Integer, Set<IAtom>> e : sphereList.entrySet()) {
//				if (e.getKey() != 0) {
//					visited.addAll(e.getValue());
//				}
//			}
			sphere.add(root);
			rankedAtomsInReactionCentre(new HashSet<IAtom>(sphere), new ArrayList<IAtom>(), 
					new ArrayList<IAtom>(sphereList.get(0)), 0, mol.getAtomCount());

			sphere = new ArrayList<IAtom>();
			sphere.addAll(sphereList.get(0));
			
			rankedAtomsOutsideReactionCentre(new HashSet<IAtom>(sphere), new ArrayList<IAtom>(), 0, 
					mol.getAtomCount()-sphere.size());
			
			Map<String,String> result = generate(sphereList, numberOfRepetitions);
			String newReactionCode = reactionCodeMapToStringOneLine(result);
//			System.out.println("newMCode " +newMCode);
			if (reactionCode == null) {
				reactionCode = newReactionCode;
				reactionCodeMap = result;
			}
			else {
				if (newReactionCode.compareTo(reactionCode) < 0) {
					reactionCode = newReactionCode;
					reactionCodeMap = result;
				}
			}
		}
		return reactionCodeMap;
	}


	/**
	 * @param sphereList
	 * @param numberOfRepetitions
	 * @return
	 */
	private Map<String,String> generate(Map<Integer, Set<IAtom>> sphereList, Map<String, Integer> numberOfRepetitions) {
		Map<Integer,String> tempResult = new TreeMap<Integer,String>();
		int leavingIndicator = sphereList.size()*2;
		
		Map<IAtom, Integer> positionsOfStayingAtoms = new HashMap<IAtom, Integer>();
		Map<IAtom, Integer> positionsOfLeavingAtoms = new HashMap<IAtom, Integer>();
		
		for (Entry<Integer, Set<IAtom>> e : sphereList.entrySet()) {
			StringBuilder stayingCodes = new StringBuilder();
			StringBuilder leavingCodes = new StringBuilder();
			String[] stayingComplementCodes = new String[3];
			String[] leavingComplementCodes = new String[3];
			
			//initialize
			stayingComplementCodes[0] = "";
			stayingComplementCodes[1] = "";
			stayingComplementCodes[2] = "";
			leavingComplementCodes[0] = "";
			leavingComplementCodes[1] = "";
			leavingComplementCodes[2] = "";
			
			int depth = e.getKey();
			List<IAtom> atoms = new ArrayList<IAtom>(e.getValue());
//			Collections.sort(atoms, new CompareByRank());
			Collections.sort(atoms, new CompareByPriority());
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
					positionsOfStayingAtoms.put(atom, positionsOfStayingAtoms.size());
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
						if (atomComplement.length() > 0) {
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
						//increment with the number of bonds (corresponds to the position of the next atom)
						cptLeaving = cptLeaving + cptBond;
					}
					else {
						stayingCodes.append(code.toString());
						
						//add complementary information (stereo and charge)
						String charges = stayingComplementCodes[0];
						String stereo = stayingComplementCodes[1];
						String isotopes = stayingComplementCodes[2];
						if (atomComplement.length() > 0) {
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
					if (atomComplement.length() > 0) {
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
	 * @param sphereList
	 * @param numberOfRepetitions
	 * @return
	 */
	/*
	private Map<String,String> generateOld(Map<Integer, Set<IAtom>> sphereList, Map<String, Integer> numberOfRepetitions) {

		Map<Integer,String> tempResult = new TreeMap<Integer,String>();
		int leavingIndicator = sphereList.size()*2;
		
		Map<IAtom, Integer> positionsOfStayingAtoms = new HashMap<IAtom, Integer>();
		Map<IAtom, Integer> positionsOfLeavingAtoms = new HashMap<IAtom, Integer>();
		
		for (Entry<Integer, Set<IAtom>> e : sphereList.entrySet()) {
			StringBuilder stayingCodes = new StringBuilder();
			StringBuilder leavingCodes = new StringBuilder();
			String[] stayingComplementCodes = new String[3];
			String[] leavingComplementCodes = new String[3];
			
			//initialize
			stayingComplementCodes[0] = "";
			stayingComplementCodes[1] = "";
			stayingComplementCodes[2] = "";
			leavingComplementCodes[0] = "";
			leavingComplementCodes[1] = "";
			leavingComplementCodes[2] = "";
			
			int depth = e.getKey();
			List<IAtom> atoms = new ArrayList<IAtom>(e.getValue());
//			Collections.sort(atoms, new CompareByRank());
			Collections.sort(atoms, new CompareByPriority());
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
					positionsOfStayingAtoms.put(atom, positionsOfStayingAtoms.size());
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
						if (atomComplement.length() > 0) {
							if (atomComplement.contains("#")) {
								int indexSymbol = atomComplement.indexOf("#");
								isotopes += String.format("%02d", cptStaying) + "-" + atomComplement.substring(indexSymbol+1, indexSymbol+3) + ";";
								atomComplement = atomComplement.substring(0, indexSymbol);
							}
							if (charge == true && stereochemistry == true && atomComplement.length() > 0) {
								if (atomComplement.contains("@")) {
									if (atomComplement.indexOf("@") > 0) {
										charges += String.format("%02d", cptLeaving) + "-" + atomComplement.substring(0, 2) + ";";
										stereo += String.format("%02d", cptLeaving) + "-" + atomComplement.substring(3) + ";";
									}
									else {
										stereo += String.format("%02d", cptLeaving) + "-" + atomComplement.substring(1) + ";";
									}
								}
								else {
									charges += String.format("%02d", cptLeaving) + "-" + atomComplement.substring(0) + ";";

								}
							}
							else if (charge == true && stereochemistry == false && atomComplement.length() > 0) {
								charges += String.format("%02d", cptLeaving) + "-" + atomComplement + ";";
							}
							else if (charge == false && stereochemistry == true && atomComplement.length() > 0) {
								stereo += String.format("%02d", cptLeaving) + "-" + atomComplement + ";";
							}
						}
						//increment by 1 because 1 atom is added
						cptLeaving++;
						for (Entry<Integer,String> e2 : bondComplements.entrySet()) {
							int index = e2.getKey() + cptLeaving;
							stereo += String.format("%02d", index) + "-" + e2.getValue() + ";";
						}
						leavingComplementCodes[0] = charges;
						leavingComplementCodes[1] = stereo;
						leavingComplementCodes[2] = isotopes;
						//increment with the number of bonds (corresponds to the position of the next atom)
						cptLeaving = cptLeaving + cptBond;
					}
					else {
						stayingCodes.append(code.toString());
						
						//add complementary information (stereo and charge)
						String charges = stayingComplementCodes[0];
						String stereo = stayingComplementCodes[1];
						String isotopes = stayingComplementCodes[2];
						if (atomComplement.length() > 0) {
							if (atomComplement.contains("#")) {
								int indexSymbol = atomComplement.indexOf("#");
								isotopes += String.format("%02d", cptStaying) + "-" + atomComplement.substring(indexSymbol+1, indexSymbol+3) + ";";
								atomComplement = atomComplement.substring(0, indexSymbol);
							}
							if (charge == true && stereochemistry == true && atomComplement.length() > 0) {
								if (atomComplement.contains("@")) {
									if (atomComplement.indexOf("@") > 0) {
										charges += String.format("%02d", cptStaying) + "-" + atomComplement.substring(0, 2) + ";";
										stereo += String.format("%02d", cptStaying) + "-" + atomComplement.substring(3) + ";";
									}
									else {
										stereo += String.format("%02d", cptStaying) + "-" + atomComplement.substring(1) + ";";
									}
								}
								else {
									charges +=String.format("%02d",  cptStaying) + "-" + atomComplement.substring(0) + ";";

								}
							}
							else if (charge == true && stereochemistry == false && atomComplement.length() > 0) {
								charges += String.format("%02d", cptStaying) + "-" + atomComplement + ";";
							}
							else if  (charge == false && stereochemistry == true && atomComplement.length() > 0) {
								stereo += String.format("%02d", cptStaying) + "-" + atomComplement + ";";
							}
						}
						//increment by 1 because 1 atom is added
						cptStaying++;
						for (Entry<Integer,String> e2 : bondComplements.entrySet()) {
							int index = e2.getKey() + cptStaying;
							stereo += String.format("%02d", index) + "-" + e2.getValue() + ";";
						}
						stayingComplementCodes[0] = charges;
						stayingComplementCodes[1] = stereo;
						stayingComplementCodes[2] = isotopes;
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
					if (atomComplement.length() > 0) {
						if (atomComplement.contains("#")) {
							int indexSymbol = atomComplement.indexOf("#");
							isotopes += String.format("%02d", cptStaying) + "-" + atomComplement.substring(indexSymbol+1, indexSymbol+3) + ";";
							atomComplement = atomComplement.substring(0, indexSymbol);
						}
						if (charge == true && stereochemistry == true && atomComplement.length() > 0) {
							if (atomComplement.contains("@")) {
								if (atomComplement.indexOf("@") > 0) {
									charges +=String.format("%02d",  cptStaying) + "-" + atomComplement.substring(0, 2) + ";";
									stereo += String.format("%02d", cptStaying) + "-" + atomComplement.substring(3) + ";";
								}
								else {
									stereo += String.format("%02d", cptStaying) + "-" + atomComplement.substring(3) + ";";
								}
							}
							else {
								charges += String.format("%02d", cptStaying) + "-" + atomComplement.substring(0) + ";";

							}
						}
						else if (charge == true && stereochemistry == false && atomComplement.length() > 0) {
							charges += String.format("%02d", cptStaying) + "-" + atomComplement + ";";
						}
						else if  (charge == false && stereochemistry == true && atomComplement.length() > 0) {
							stereo += String.format("%02d", cptStaying) + "-" + atomComplement.replace("@","") + ";";
						}
						
					}
					//increment by 1 because 1 atom is added
					cptStaying++;
					for (Entry<Integer,String> e2 : bondComplements.entrySet()) {
						int index = e2.getKey() + cptStaying;
						stereo += String.format("%02d", index) + "-" + e2.getValue() + ";";
					}
					stayingComplementCodes[0] = charges;
					stayingComplementCodes[1] = stereo;
					stayingComplementCodes[2] = isotopes;
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
	*/

	/**
	 * @param acSet
	 * @return
	 */
	private HashMap<String, IAtom> indexAtomIAtomContainerSet (IAtomContainerSet acSet) {
		HashMap<String, IAtom> result = new HashMap<String, IAtom>();

		for (int i = 0; i < acSet.getAtomContainerCount(); i ++) {
			IAtomContainer ac = acSet.getAtomContainer(i);
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
				chargeInReactant = encodeChargeAndIsotope(atomInReactant.getFormalCharge());
				chargeInProduct = encodeChargeAndIsotope(atomInProduct.getFormalCharge());
			}
			else {
				chargeInReactant = "0";
				chargeInProduct = encodeChargeAndIsotope(atomInProduct.getFormalCharge());
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
		
		String isotopeInReactant = "0";
		String isotopeInProduct = "0";
		if (atomInProduct.getMassNumber() != null && !atomInProduct.getSymbol().equals("H")) {
			int frequentMass = ElementCalculation.calculateMass(atomInProduct.getSymbol());
			if (frequentMass != atomInProduct.getMassNumber()) {
				int delta = atomInProduct.getMassNumber() - frequentMass;
				isotopeInProduct = encodeChargeAndIsotope(delta);
			}
		}
		if (atomInReactant.getMassNumber() != null && !atomInReactant.getSymbol().equals("H")) {
			int frequentMass = ElementCalculation.calculateMass(atomInReactant.getSymbol());
			if (frequentMass != atomInReactant.getMassNumber()) {
				int delta = atomInReactant.getMassNumber() - frequentMass;
				isotopeInReactant = encodeChargeAndIsotope(delta);
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
			if ((int) inProduct.getProperty(additionalConstants.BOND_CHANGE_INFORMATION) == additionalConstants.BOND_FORMED) {
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
	
//	/**
//	 * @param bond
//	 * @return
//	 */
//	private String encodeBondStereoChange(IBond bond) {
//		if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE) == null) {
//			return "0";
//		}
//		else {
//			if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.NONE_TO_UP)) {
//				return "A";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.NONE_TO_DOWN)) {
//				return "B";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.NONE_TO_UP_INVERTED)) {
//				return "C";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.NONE_TO_DOWN_INVERTED)) {
//				return "D";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.NONE_TO_Z)) {
//				return "E";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.NONE_TO_E)) {
//				return "F";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.UP_TO_NONE)) {
//				return "G";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.DOWN_TO_NONE)) {
//				return "H";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.UP_INVERTED_TO_NONE)) {
//				return "I";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.DOWN_INVERTED_TO_NONE)) {
//				return "J";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.Z_TO_NONE)) {
//				return "K";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.E_TO_NONE)) {
//				return "L";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.DOWN_TO_UP)) {
//				return "Z";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.UP_TO_DOWN)) {
//				return "W";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.DOWN_INVERTED_TO_UP_INVERTED)) {
//				return "Y";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.UP_INVERTED_TO_DOWN_INVERTED)) {
//				return "X";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.E_TO_Z)) {
//				return "Q";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.Z_TO_E)) {
//				return "R";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.NONE_TO_UP)) {
//				return "M";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.NONE_TO_DOWN)) {
//				return "N";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.NONE_TO_UP_INVERTED)) {
//				return "O";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.NONE_TO_DOWN_INVERTED)) {
//				return "P";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.NONE_TO_Z)) {
//				return "9";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.NONE_TO_E)) {
//				return "8";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.DOWN__TO_UP_INVERTED)) {
//				return "V";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.UP_TO_DOWN_INVERTED)) {
//				return "U";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.DOWN_INVERTED_TO_UP)) {
//				return "T";
//			}
//			else if (bond.getProperty(additionalConstants.BOND_STEREO_CHANGE).equals(additionalConstants.UP_INVERTED_TO_DOWN)) {
//				return "S";
//			}
//		}
//		return "0";
//	}


	/**
	 * @param chargeOrMass
	 * @return
	 */
	private String encodeChargeAndIsotope(int chargeOrMass) {
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
			if (atom.getProperty("depth") == null) {
				atom.setProperty("depth", depth);
			}
			List<IBond> bonds = atomContainer.getConnectedBondsList(atom);
			if (atom.getProperty("code") == null) {
				String code = encodeAtom(atom, bonds, indexAtomReactants);
				atom.setProperty("code", code+"_"+cpt);
			}
			for (IBond bond : bonds) {
				if (bond.getProperty("code") == null) {
//					bond.setProperty("code", encodeBond(bond)+"_"+depth);
					bond.setProperty("code", encodeBond(bond, indexBondReactants)+"_"+depth);
				}
				nextAtom = bond.getOther(atom);
				if (!visited.contains(nextAtom)) {
					if (!sphere.contains(nextAtom)) {
						newSphere.add(nextAtom);
					}
				}
			}
		}
		if (newSphere.size() > 0) {
			if (depth == 0) {
				sphereList.put(0, newSphere);
			}
			sphereList.put(depth, newSphere);
			//System.out.println(sphereList);
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

//	/**
//	 * @param hybridization
//	 * @return
//	 */
//	private IAtomType.Hybridization hybridizationDecoding(String hybridization) {
//		if (Integer.parseInt(hybridization,16) == 1) return IAtomType.Hybridization.PLANAR3;
//		if (Integer.parseInt(hybridization,16) == 2) return IAtomType.Hybridization.S;
//		if (Integer.parseInt(hybridization,16) == 3) return IAtomType.Hybridization.SP1;
//		if (Integer.parseInt(hybridization,16) == 4) return IAtomType.Hybridization.SP2;
//		if (Integer.parseInt(hybridization,16) == 5) return IAtomType.Hybridization.SP3;
//		if (Integer.parseInt(hybridization,16) == 6) return IAtomType.Hybridization.SP3D1;
//		if (Integer.parseInt(hybridization,16) == 7) return IAtomType.Hybridization.SP3D2;
//		if (Integer.parseInt(hybridization,16) == 8) return IAtomType.Hybridization.SP3D3;
//		if (Integer.parseInt(hybridization,16) == 9) return IAtomType.Hybridization.SP3D4;
//		if (Integer.parseInt(hybridization,16) == 10) return IAtomType.Hybridization.SP3D5;
//		return null;
//	}

//	/**
//	 * @param acSet
//	 * @return
//	 */
//	private HashMap<String, IAtomType.Hybridization> indexIDHybridization (IAtomContainerSet acSet) {
//		HashMap<String, IAtomType.Hybridization> result = new HashMap<String, IAtomType.Hybridization>();
//
//		for (int i = 0; i < acSet.getAtomContainerCount(); i ++) {
//			IAtomContainer ac = acSet.getAtomContainer(i);
//			for (int j = 0; j < ac.getAtomCount(); j++) {
//				result.put(ac.getAtom(j).getID(), ac.getAtom(j).getHybridization());
//			}
//		}
//		return result;
//	}

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
	 * @param bond
	 * @return
	 */
	private int getRCScore(List<IBond> bonds) {
		int centerStatusValue = 99;
		for (IBond bond : bonds) {
			int score = 0;
			if (bond.getProperty(BOND_CHANGE_INFORMATION) != null) {
				score += 1;
				if ((int) bond.getProperty(BOND_CHANGE_INFORMATION) == additionalConstants.BOND_CLEAVED) 
					score += 6;
				if ((int) bond.getProperty(BOND_CHANGE_INFORMATION) == additionalConstants.BOND_FORMED) 
					score += 8;
				if ((int) bond.getProperty(BOND_CHANGE_INFORMATION) == additionalConstants.BOND_ORDER) 
					score += 4;
			}
			else {
				score += 0;
			}
			if (score > 0 && score < centerStatusValue) {
				centerStatusValue = score;
			}
		}
		
		return (centerStatusValue != 99) ? centerStatusValue : 0;
	}

	/**
	 * @param atomInRC
	 * @return
	 */
	private List<IAtom> getRootAtoms(List<IAtom> atomInRC) {
		List<IAtom> root = new ArrayList<IAtom>();
		if (atomInRC.size() < 2) {
			root.add(atomInRC.get(0));
			return root;
		}
		else {
			if (!atomInRC.get(0).getProperty("code").equals(atomInRC.get(1).getProperty("code"))) {
				root.add(atomInRC.get(0));
				return root;
			}
			//compare
			else {
				List<IAtom> solutions = new ArrayList<IAtom>();
				Map<Integer, Set<IAtom>> conflicts = new HashMap<Integer, Set<IAtom>>();
				Map<Integer,Integer> colors = new HashMap<Integer,Integer>();
				for (int i = 0; i < atomInRC.size(); i++) {
					if (i == 0) {
						solutions.add(atomInRC.get(i));
						Set<IAtom> set = new HashSet<IAtom>();
						set.add(atomInRC.get(i));
						conflicts.put(i, set);
						colors.put(i, 1);
					}
					else if (atomInRC.get(0).getProperty("code").equals(atomInRC.get(i).getProperty("code"))){
						solutions.add(atomInRC.get(i));
						Set<IAtom> set = new HashSet<IAtom>();
						set.add(atomInRC.get(i));
						conflicts.put(i, set);
						colors.put(i, 1);
					}
					else{
						break;
					}
				}
				//all path are similar ex c1ccccc1
				if (solutions.size() == atomInRC.size()) {
					root.add(atomInRC.get(0));
					return root;
				}
				//Explore neighbors to get best solution
				else {
					//C1CCN2CCCCN2C1
//					Map<IAtom,Integer> colors = initColors(solutions);
//					attributeColorsByNextNeighbors(solutions, colors, new HashMap<IAtom,List<IAtom>>(), new HashMap<IAtom,List<IAtom>>());
					attributeColorsByNextNeighbors(conflicts, colors, 
							new HashMap<Integer, String>(), new HashMap<Integer,List<IAtom>>()) ;
//					System.out.println(colors);
//					return indexColorsToList(colors).get(0);
					root = makeListOfAtoms(getIndexesOfAtomsWithTheBestColor(colors), atomInRC);
					return root;
				}
			}
		}
	}
	
	/**
	 * @param indexes
	 * @param atomInRC
	 * @return
	 */
	private List<IAtom> makeListOfAtoms(List<Integer> indexes, List<IAtom> atomInRC) {
		List<IAtom> res = new ArrayList<IAtom>();
		for (int index : indexes) {
			res.add(atomInRC.get(index));
		}
		return res;
	}
	
	/**
	 * @param colors
	 * @return
	 */
	private List<Integer> getIndexesOfAtomsWithTheBestColor(Map<Integer,Integer> colors) {
		List<Integer> result = new ArrayList<Integer>();
		
		List<Integer> allColors = new ArrayList<Integer>(colors.values());
		Collections.sort(allColors);
		
		int bestColor = allColors.get(0);
		
		for (Entry<Integer,Integer> e : colors.entrySet()) {
			int color = e.getValue();
			
			if (color == bestColor) {
				result.add(e.getKey());
			}
		}

		return result;
	}
	
	/**
	 * @param conflicts
	 * @param colorOfAtoms
	 * @param paths
	 * @param visited
	 */
	private void attributeColorsByNextNeighbors(Map<Integer, Set<IAtom>> conflicts, Map<Integer,Integer> colorOfAtoms, 
			Map<Integer, String> paths, Map<Integer,List<IAtom>> visited) {
		
		if (conflicts.isEmpty()) {
			return;
		}
		else {
			//update path code
			List<Integer> keysToRemoveInConflicts = new ArrayList<Integer>(); 
			for (Entry<Integer, Set<IAtom>> e : conflicts.entrySet()) {
				int key = e.getKey();
				Set<IAtom> newSphere = new HashSet<IAtom>();
				Set<IAtom> sphere = e.getValue();
				List<IAtom> visitedAtoms;
				if (visited.get(key) != null) {
					visitedAtoms = visited.get(key);
				}
				else {
					visitedAtoms = new ArrayList<IAtom>();
				}
				
				List<String> pathCodes = new ArrayList<String>();
//				System.out.println(sphere.size() + " " + conflicts.size());
				visitedAtoms.addAll(sphere);
				for (IAtom atom : sphere) {
					List<IBond> bonds = mol.getConnectedBondsList(atom);
					for (IBond bond : bonds) {
						IAtom other = bond.getOther(atom);
						if (!visitedAtoms.contains(other)) {
							newSphere.add(other);
							StringBuilder codeBuilder = new StringBuilder();
							codeBuilder.append(other.getProperty("code").toString()+ "-");
							List<IBond> bondsOfOther = mol.getConnectedBondsList(other);
							List<String> bondsCode = new ArrayList<String>();
							for (IBond bondOfOther : bondsOfOther) {
								bondsCode.add(bondOfOther.getProperty("code").toString());
							}
							Collections.sort(bondsCode, Collections.reverseOrder());
							for (String bcode : bondsCode) {
								codeBuilder.append(bcode);
							}
							codeBuilder.append("_");
							pathCodes.add(codeBuilder.toString());
						}
					}
				}
				if (!newSphere.isEmpty()) {
					conflicts.put(key, newSphere);
				}
				else {
					keysToRemoveInConflicts.add(key);
				}
				visited.put(key, visitedAtoms);
				
				Collections.sort(pathCodes, Collections.reverseOrder());
				StringBuilder path = new StringBuilder();
				for (String s : pathCodes) {
					path.append(s);
				}
				String newPath; 
				if (paths.containsKey(key)) {
					newPath = paths.get(key) + path.toString();
				}
				else {
					newPath = path.toString();
				}
//				System.out.println(key + " " +newPath);
				paths.put(key, newPath);
			}
			
			//update colors
			Map<Integer,Integer> colorOfAtomsCopy = new HashMap<Integer,Integer>(colorOfAtoms);
			for (int i : paths.keySet()) {
				String path1 = paths.get(i);
//				System.out.println("p1 "+ path1);
				if (!conflicts.containsKey(i)) {
					continue;
				}
				for (int j : paths.keySet()) {
					if (i == j) {
						continue;
					}
					if (!conflicts.containsKey(j)) {
						continue;
					}
					String path2 = paths.get(j);
//					System.out.println("p1 "+ i + " " + path1);
//					System.out.println("p2 "+ j + " " + path2);
					if (path1.compareTo(path2) > 0 && colorOfAtomsCopy.get(i) == colorOfAtomsCopy.get(j)) {
						int color = colorOfAtoms.get(j) + 1;
						colorOfAtoms.put(j, color);
					}
				}
			}
			
			//find index of atom(s) with an unique color and remove it of the next iteration
			Map<Integer,List<Integer>> sameColorCounter = new HashMap<Integer,List<Integer>>();
			for (int i = 0; i < colorOfAtoms.size(); i++) {
				int color = colorOfAtoms.get(i);
				if (!sameColorCounter.containsKey(color)) {
					List<Integer> indexes = new ArrayList<Integer>();
					indexes.add(i);
					sameColorCounter.put(color, indexes);
				}
				else {
					List<Integer> indexes = sameColorCounter.get(color);
					indexes.add(i);
					sameColorCounter.put(color, indexes);
				}
			}
			for (List<Integer> list : sameColorCounter.values()) {
				if (list.size() == 1) {
					int key = list.get(0);
//					System.out.println(key + " " + paths.get(key));
					conflicts.remove(key);
					paths.remove(key);
					visited.remove(key);
				}
			}
			
			for (int k : keysToRemoveInConflicts) {
				conflicts.remove(k);
			}
//			System.out.println("colors " +colorOfAtoms);
			attributeColorsByNextNeighbors(conflicts, colorOfAtoms, paths, visited);
		}
	}
	
	/**
	 * 
	 */
	private void resetRankingProperties() {
		for (IAtom atom : mol.atoms()) {
			atom.removeProperty("priority");
			atom.removeProperty("rank");
			atom.removeProperty("pos");
		}
	}
	
	/**
	 * @param atom
	 */
	private void resetRankingProperties(IAtom atom) {			
		atom.removeProperty("priority");
		atom.removeProperty("rank");
		atom.removeProperty("pos");
	}
	
	/**
	 * Attribute a priority of each atom in the reaction centre in order to be later ranked. 
	 * NB: List to process contains the atom in the reaction center and the algorithm continue until this list is not null
	 * (ie all atoms in the RC have been processed). Indeed, 2 independent reaction center in the pseudo molecule may occur,
	 * in that case, the algorithm need to classify in function of the determined root.
	 * @param sphere
	 * @param visited
	 * @param toProcess
	 * @param depth
	 * @param pos
	 * @return
	 */
	private void rankedAtomsInReactionCentre(Set<IAtom> sphere, List<IAtom> visited, 
			List<IAtom> toProcess, int depth, int pos) {
		Set<IAtom> newSphere = new HashSet<IAtom>();
		Map<IAtom, IAtom> prevAdj = new HashMap<IAtom, IAtom>();
		Map<String, List<IAtom>> codesOfAtoms = new TreeMap<String, List<IAtom>>(Collections.reverseOrder()); 
		
		if (depth == 0) {
			IAtom first = sphere.iterator().next();
			first.setProperty("priority", first.getProperty("code") +  "_" + pos);
			first.setProperty("rank",depth + "_" + 0);
			first.setProperty("pos", pos);
			depth++;
			pos--;
		}

		visited.addAll(sphere);
		toProcess.removeAll(sphere);
		for (IAtom atom : sphere) {
			List<IBond> bonds = mol.getConnectedBondsList(atom);
			for (IBond bond : bonds) {
				IAtom other = bond.getOther(atom);
				if (!visited.contains(other)) {
					String code = depth + "_" + other.getProperty("code").toString() + "-";
					List<IBond> bondsOfOther = mol.getConnectedBondsList(other);
					List<String> bondsCode = new ArrayList<String>();
					int posCounter = 0;
					for (IBond bondOfOther : bondsOfOther) {
						bondsCode.add(bondOfOther.getProperty("code").toString());
						IAtom other2 = bondOfOther.getOther(other);
						if (other2.getProperty("pos") != null) {
							int cPos = other2.getProperty("pos");
							posCounter += cPos;
						}
					}
					code += ":" + posCounter;
					Collections.sort(bondsCode, Collections.reverseOrder());
					for (String bcode : bondsCode) {
						code += bcode;
					}
					
					if (!codesOfAtoms.containsKey(code)) {
						List<IAtom> list = new ArrayList<IAtom>();
						list.add(other);
						codesOfAtoms.put(code, list);
					}
					else {
						List<IAtom> list = codesOfAtoms.get(code);
						if (!list.contains(other)) {
							list.add(other);
							codesOfAtoms.put(code, list);
						}
					}
					if (!visited.contains(other)) {
						newSphere.add(other);
					}
				}
				if (prevAdj.containsKey(other)) {
					IAtom prev = prevAdj.get(other);
					if (atom.getProperty("rank").toString().compareTo(prev.getProperty("rank").toString()) > 0) {
						prevAdj.put(other, prev);
					}
				}
				else {
					prevAdj.put(other, atom);
				}
			}
		}
		
		int cpt = 0;
		List<IAtom> processed = new ArrayList<IAtom>();
//		System.out.println(codesOfAtoms.keySet());
		for (Entry<String, List<IAtom>> e : codesOfAtoms.entrySet()) {
//			String code = e.getKey();
			List<IAtom> atoms = e.getValue();
//			System.out.println("codesOfAtoms " + code + " atomSize " + atoms.size());
//			System.out.println("cpt " + cpt + " codesOfAtoms " + code + " atomSize " + atoms.size());
			if (atoms.size() == 1 && !processed.contains(atoms.get(0))) {
//				code += cpt;
				IAtom atom = atoms.get(0);
				atom.setProperty("priority", atom.getProperty("code") +  "_" + pos);
				atom.setProperty("rank", depth + "_" + cpt);
				atom.setProperty("pos", pos);
//				System.out.println("depth " + depth + " " + code + " " + atoms.get(0).getProperty("rank") + " " +
//						atoms.get(0).getProperty("uID"));
				processed.add(atoms.get(0));
				pos--;
				cpt++;
			}
			else {
				//classify
//				if (atoms.get(0).getProperty("rank") != null) {
//					//comparator according to rank
//				}
//				else {
					//solve using previous
					Map<Integer, Set<IAtom>> conflicts = new HashMap<Integer, Set<IAtom>>();
					Map<Integer,Integer> colors = new HashMap<Integer,Integer>();
					for (int i = 0; i < atoms.size(); i++) {
						Set<IAtom> set = new HashSet<IAtom>();
						set.add(atoms.get(i));
						conflicts.put(i, set);
						colors.put(i, 0);
//						System.out.println(atoms.get(i).getProperties());
					}
					
					//update colors
					Map<Integer,Integer> colorOfAtomsCopy = new HashMap<Integer,Integer>(colors);
					for (int i = 0; i < atoms.size(); i++) {
//						System.out.println(prevAdj.get(atoms.get(i)).getProperties());
						String rank1 = prevAdj.get(atoms.get(i)).getProperty("rank");
						for (int j = 0; j < atoms.size(); j++) {
							if (i == j) {
								continue;
							}
							String rank2 = prevAdj.get(atoms.get(j)).getProperty("rank");
							if (rank1.compareTo(rank2) < 0 && colorOfAtomsCopy.get(i) == colorOfAtomsCopy.get(j)) {
								int color = colors.get(j) + 1;
								colors.put(j, color);
							}
						}
					}
					
					//find index of atom(s) with an unique color and remove it of the next iteration
					Map<Integer,List<Integer>> sameColorCounter = new HashMap<Integer,List<Integer>>();
					for (int i = 0; i < colors.size(); i++) {
						int color = colors.get(i);
						if (!sameColorCounter.containsKey(color)) {
							List<Integer> indexes = new ArrayList<Integer>();
							indexes.add(i);
							sameColorCounter.put(color, indexes);
						}
						else {
							List<Integer> indexes = sameColorCounter.get(color);
							indexes.add(i);
							sameColorCounter.put(color, indexes);
						}
					}
					int cpt2 = 0;
					for (List<Integer> list : sameColorCounter.values()) {
						if (list.size() == 1) {
							int key = list.get(0);
							conflicts.remove(key);
							IAtom atom = atoms.get(key);
							int score = cpt + colors.get(key);
							atom.setProperty("priority", atom.getProperty("code") +  "_" + pos);
							atom.setProperty("rank", depth + "_" + score);
							atom.setProperty("pos", pos);
//							System.out.println("depth " + depth + " " + code + " " +atom.getProperty("rank") + 
//									 " " +atom.getProperty("uID"));
							processed.add(atom);
							pos--;
							cpt2++;
						}
					}

					//solve using next
					if (conflicts.size() > 0) {
//						System.out.println("depth " + depth + " " +conflicts);
//						for (int l = 0; l < conflicts.size(); l++) {
//							System.out.println(conflicts.get(l).iterator().next().getProperties() + " " + conflicts.get(l).iterator().next());
//						}
						//Explore neighbors to get best solution
						attributeColorsByNextNeighbors(conflicts, colors, 
									new HashMap<Integer, String>(), new HashMap<Integer,List<IAtom>>());
					}
//					System.out.println(colors);
					
					//ascending sort by color 
					Map<Integer,Integer> sortedAtomIndexByColors = sortByValue(colors);
					
//					System.out.println("sort " + sortedAtomIndexByColors + " " + cpt);
					
					for (Entry<Integer,Integer> e2 : sortedAtomIndexByColors.entrySet()) {
						IAtom atom = atoms.get(e2.getKey());
						int color = e2.getValue();
						if (!processed.contains(atom)) {
							int score = cpt + color;
//							System.out.println("score " + score + " cpt " + cpt + " color "+ color);
							atom.setProperty("priority", atom.getProperty("code") +  "_" + pos);
							atom.setProperty("rank", depth + "_" + score);
							atom.setProperty("pos", pos);
//							System.out.println("depth " + depth + " " + code + " " +atom.getProperty("rank") + 
//									 " " +atom.getProperty("uID"));							
							processed.add(atom);
							pos--;
							cpt2++;
						}
					}
					cpt += cpt2;
//				}
				
				
				//1 look in first atom if rank exist
				//2 if not compare using best prev
				//3 if not compare using next neighboors by making all path and then set rank
					//map<IAtom[], code> ne pas oublier de mettre copy de visited pour exclure les atomes deja visites
				//attirbuer color type depth+codeBOnd+mcode
			}
		}
		
//		for (IAtom atom : newSphere) {
//			System.out.println("depth " + depth + " code " + atom.getProperty("code") + " rank " + atom.getProperty("rank") + " uID " +
//					atom.getProperty("uID"));
//		}
		
		if (newSphere.isEmpty() || toProcess.isEmpty()) {
			//remove property of non reaction centre atoms
			for (IAtom atom : visited) {
				if (atom.getProperty(additionalConstants.REACTION_CENTER) == null) {
					resetRankingProperties(atom);
				}
			}
			return;
		}
		else {
			rankedAtomsInReactionCentre(newSphere, visited, toProcess, depth+1, pos);
		}
	}

	/**
	 * @param sphere
	 * @param visited
	 * @param depth
	 * @param pos
	 */
	private void rankedAtomsOutsideReactionCentre(Set<IAtom> sphere, List<IAtom> visited, int depth, int pos) {
		Set<IAtom> newSphere = new HashSet<IAtom>();
		Map<IAtom, IAtom> prevAdj = new HashMap<IAtom, IAtom>();
		Map<String, List<IAtom>> codesOfAtoms = new TreeMap<String, List<IAtom>>(Collections.reverseOrder()); 
		
		visited.addAll(sphere);
		for (IAtom atom : sphere) {
			List<IBond> bonds = mol.getConnectedBondsList(atom);
			for (IBond bond : bonds) {
				IAtom other = bond.getOther(atom);
				if (!visited.contains(other)) {
					String code = depth + "_" + other.getProperty("code").toString() + "-";
					List<IBond> bondsOfOther = mol.getConnectedBondsList(other);
					List<String> bondsCode = new ArrayList<String>();
					int posCounter = 0;
					for (IBond bondOfOther : bondsOfOther) {
						bondsCode.add(bondOfOther.getProperty("code").toString());
						IAtom other2 = bondOfOther.getOther(other);
						if (other2.getProperty("pos") != null) {
							int cPos = other2.getProperty("pos");
							posCounter += cPos;
						}
					}
					code += ":" + posCounter;
					Collections.sort(bondsCode, Collections.reverseOrder());
					for (String bcode : bondsCode) {
						code += bcode;
					}
					
					if (!codesOfAtoms.containsKey(code)) {
						List<IAtom> list = new ArrayList<IAtom>();
						list.add(other);
						codesOfAtoms.put(code, list);
					}
					else {
						List<IAtom> list = codesOfAtoms.get(code);
						if (!list.contains(other)) {
							list.add(other);
							codesOfAtoms.put(code, list);
						}
					}
					newSphere.add(other);
				}
				if (prevAdj.containsKey(other)) {
					IAtom prev = prevAdj.get(other);
					if (atom.getProperty("rank").toString().compareTo(prev.getProperty("rank").toString()) > 0) {
						prevAdj.put(other, prev);
					}
				}
				else {
					prevAdj.put(other, atom);
				}
			}
		}
		
		int cpt = 0;
		List<IAtom> processed = new ArrayList<IAtom>();
//		System.out.println(codesOfAtoms.keySet());
		for (Entry<String, List<IAtom>> e : codesOfAtoms.entrySet()) {
//			String code = e.getKey();
			List<IAtom> atoms = e.getValue();
//			System.out.println("codesOfAtoms " + code + " atomSize " + atoms.size());
//			System.out.println("cpt " + cpt + " codesOfAtoms " + code + " atomSize " + atoms.size());
			if (atoms.size() == 1 && !processed.contains(atoms.get(0))) {
//				code += cpt;
				IAtom atom = atoms.get(0);
				atom.setProperty("priority", atom.getProperty("code") +  "_" + pos);
				atom.setProperty("rank", depth + "_" + cpt);
				atom.setProperty("pos", pos);
//				System.out.println("depth " + depth + " " + code + " " + atoms.get(0).getProperty("rank") + " " +
//						atoms.get(0).getProperty("uID"));
				processed.add(atoms.get(0));
				pos--;
				cpt++;
			}
			else {
				//classify
//				if (atoms.get(0).getProperty("rank") != null) {
//					//comparator according to rank
//				}
//				else {
					//solve using previous
					Map<Integer, Set<IAtom>> conflicts = new HashMap<Integer, Set<IAtom>>();
					Map<Integer,Integer> colors = new HashMap<Integer,Integer>();
					for (int i = 0; i < atoms.size(); i++) {
						Set<IAtom> set = new HashSet<IAtom>();
						set.add(atoms.get(i));
						conflicts.put(i, set);
						colors.put(i, 0);
//						System.out.println(atoms.get(i).getProperties());
					}
					
					//update colors
					Map<Integer,Integer> colorOfAtomsCopy = new HashMap<Integer,Integer>(colors);
					for (int i = 0; i < atoms.size(); i++) {
//						System.out.println(prevAdj.get(atoms.get(i)).getProperties());
						String rank1 = prevAdj.get(atoms.get(i)).getProperty("rank");
						for (int j = 0; j < atoms.size(); j++) {
							if (i == j) {
								continue;
							}
							String rank2 = prevAdj.get(atoms.get(j)).getProperty("rank");
							if (rank1.compareTo(rank2) < 0 && colorOfAtomsCopy.get(i) == colorOfAtomsCopy.get(j)) {
								int color = colors.get(j) + 1;
								colors.put(j, color);
							}
						}
					}
					
					//find index of atom(s) with an unique color and remove it of the next iteration
					Map<Integer,List<Integer>> sameColorCounter = new HashMap<Integer,List<Integer>>();
					for (int i = 0; i < colors.size(); i++) {
						int color = colors.get(i);
						if (!sameColorCounter.containsKey(color)) {
							List<Integer> indexes = new ArrayList<Integer>();
							indexes.add(i);
							sameColorCounter.put(color, indexes);
						}
						else {
							List<Integer> indexes = sameColorCounter.get(color);
							indexes.add(i);
							sameColorCounter.put(color, indexes);
						}
					}
					int cpt2 = 0;
					for (List<Integer> list : sameColorCounter.values()) {
						if (list.size() == 1) {
							int key = list.get(0);
							conflicts.remove(key);
							IAtom atom = atoms.get(key);
							int score = cpt + colors.get(key);
							atom.setProperty("priority", atom.getProperty("code") +  "_" + pos);
							atom.setProperty("rank", depth + "_" + score);
							atom.setProperty("pos", pos);
//							System.out.println("depth " + depth + " " + code + " " +atom.getProperty("rank") + 
//									 " " +atom.getProperty("uID"));
							processed.add(atom);
							pos--;
							cpt2++;
						}
					}

					//solve using next
					if (conflicts.size() > 0) {
//						System.out.println("depth " + depth + " " +conflicts);
//						for (int l = 0; l < conflicts.size(); l++) {
//							System.out.println(conflicts.get(l).iterator().next().getProperties() + " " + conflicts.get(l).iterator().next());
//						}
						//Explore neighbors to get best solution
						attributeColorsByNextNeighbors(conflicts, colors, 
									new HashMap<Integer, String>(), new HashMap<Integer,List<IAtom>>());
					}
//					System.out.println(colors);
					
					//ascending sort by color 
					Map<Integer,Integer> sortedAtomIndexByColors = sortByValue(colors);
					
//					System.out.println("sort " + sortedAtomIndexByColors + " " + cpt);
					
					for (Entry<Integer,Integer> e2 : sortedAtomIndexByColors.entrySet()) {
						IAtom atom = atoms.get(e2.getKey());
						int color = e2.getValue();
						if (!processed.contains(atom)) {
							int score = cpt + color;
//							System.out.println("score " + score + " cpt " + cpt + " color "+ color);
							atom.setProperty("priority", atom.getProperty("code") +  "_" + pos);
							atom.setProperty("rank", depth + "_" + score);
							atom.setProperty("pos", pos);
//							System.out.println("depth " + depth + " " + code + " " +atom.getProperty("rank") + 
//									 " " +atom.getProperty("uID"));							
							processed.add(atom);
							pos--;
							cpt2++;
						}
					}
					cpt += cpt2;
//				}
				

			}
		}
		
//		for (IAtom atom : newSphere) {
//			System.out.println("depth " + depth + " code " + atom.getProperty("code") + " rank " + atom.getProperty("rank") + " uID " +
//					atom.getProperty("uID"));
//		}
		
		if (newSphere.isEmpty()) {
			return;
		}
		else {
			rankedAtomsOutsideReactionCentre(newSphere, visited, depth+1, pos);
		}
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
	private <K, V extends Comparable<? super V>> Map<K, V> sortByValue(Map<K, V> map) {
        List<Entry<K, V>> list = new ArrayList<>(map.entrySet());
        list.sort(Entry.comparingByValue());

        Map<K, V> result = new LinkedHashMap<>();
        for (Entry<K, V> entry : list) {
            result.put(entry.getKey(), entry.getValue());
        }

        return result;
    }
	
	/**
	 * @return
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
	private void reinitializedConstants(IAtomContainer container) {
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
}

//descending comparison res = (3,2,1)
class CompareByPriority implements Comparator<IAtom> { 
	public int compare(IAtom a1, IAtom a2) { 
		if (((String) a1.getProperty("priority")).compareTo(a2.getProperty("priority")) < 0)
			return 1;
		else 
			return -1;
	} 
}

class CompareByRank implements Comparator<IAtom>{ 
	public int compare(IAtom a1, IAtom a2) { 
		if ((a1.getProperty("rank").toString()).compareTo(a2.getProperty("rank").toString()) < 0)
			return -1;
		else 
			return 1;
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
