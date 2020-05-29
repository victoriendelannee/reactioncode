/**
 * 
 */
package com.nih.reaction;

import static com.nih.reaction.additionalConstants.BOND_CHANGE_INFORMATION;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import static java.io.File.separator;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.Bond;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality;
import org.openscience.cdk.layout.StructureDiagramGenerator2;
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
//import org.openscience.smsd.tools.ExtAtomContainerManipulator;

import com.nih.writer.MDLV2000Writer;

/**
 * @author delanneev
 *
 */
public class PseudoMolecule {

	IAtomContainer pseudoMolecule;
	final Set<IBond> bondFormedList = new HashSet<IBond>();
	final Set<IBond> bondCleavedList = new HashSet<IBond>();
	Set<IBond> bondOrderList = new HashSet<IBond>();
	Set<IBond> bondStereoList = new HashSet<IBond>();
	final Set<IAtom> reactionCenter = new HashSet<IAtom>();
	HashMap<String, Integer> atomRepetitions = new HashMap<String, Integer>();
	Map<String,LinkedHashSet<IAtom>> indexAtomsPseudoMol = new HashMap<String,LinkedHashSet<IAtom>>();
	Map<String,LinkedHashSet<IAtom>> indexAtomsPseudoMolUnique = new HashMap<String,LinkedHashSet<IAtom>>();

	Set<String> atomsInPseudoMol = new HashSet<String>();
	Set<String> bondsInPseudoMol = new HashSet<String>();
	Set<IAtom> leavingAtom = new HashSet<IAtom>();
	Set<IBond> leavingBond = new HashSet<IBond>();

	boolean repetitionToProcess = false;

	/**
	 * @param args
	 */
	public void main(String[] args) {
		// TODO Auto-generated method stub

	}

	//new version
	/**
	 * Make the pseudo-molecule
	 * @param reactants
	 * @param products
	 * @param builder
	 * @return
	 * @throws CDKException 
	 */
	public IAtomContainer makePseudoMolecule(IAtomContainerSet reactants, IAtomContainerSet products) throws CDKException {

		bondFormedList.clear();
		bondCleavedList.clear();
		reactionCenter.clear();
		atomsInPseudoMol.clear();
		bondsInPseudoMol.clear();
		atomRepetitions.clear();
		indexAtomsPseudoMolUnique.clear();	
		repetitionToProcess = false;

		pseudoMolecule = DefaultChemObjectBuilder.getInstance().newAtomContainer();

		leavingAtom = new HashSet<IAtom>();
		leavingBond = new HashSet<IBond>();
		
		//build products atom index and manage potential repeated atoms
		IAtomContainer temp = buildIndexInPseudoMolecule(products);
		if (repetitionToProcess) 
			duplicateAtomManagement(temp);

		for (IAtomContainer product : products.atomContainers()) { 
			for (IAtom a : product.atoms()) {
				addAtomInPseudoMolecule(indexAtomsPseudoMolUnique.get(a.getID()).iterator().next(), 0);
			}
			for (IBond b : product.bonds()) {
					// Parameterize and add change bond Information
					addBondInPseudoMolecule(b, 0);
			}
		}
		for (IAtomContainer reactant : reactants.atomContainers()) { 

			for (IAtom a : reactant.atoms()) {
				//pseudoMolecule.addAtom(reactant.getAtom(j));
				//indexAtomsReactants.put(reactant.getAtom(j).getID(),reactant.getAtom(j));
				if (indexAtomsPseudoMolUnique.containsKey(a.getID()) == false) {
					addAtomInPseudoMolecule(a, 1);
				}

			}
			for (IBond b : reactant.bonds()) {
				addBondInPseudoMolecule(b, 1);
			}
		}
		
		//DO NOT AROMATIZE because it could modify the right aromaticity which is already set up
		//ExtAtomContainerManipulator.aromatizeMolecule(pseudoMolecule);

		return cleanMolecule(pseudoMolecule);
	}

	/**
	 * @param atom
	 * @param leaving
	 */
	private void addAtomInPseudoMolecule(IAtom atom, int leaving) {
		if (!atomsInPseudoMol.contains(atom.getID())) {
			pseudoMolecule.addAtom(atom);
			atomsInPseudoMol.add(atom.getID());
			
			//need to add missing atom in indexAtomsPseudoMol (= atoms in reactants but no in products)
			if (indexAtomsPseudoMolUnique.get(atom.getID()) == null) {
				repetitionToProcess = true;
				LinkedHashSet<IAtom> t = new LinkedHashSet<IAtom>();
				t.add(atom);
				indexAtomsPseudoMolUnique.put(atom.getID(), t);
				indexAtomsPseudoMol.put(atom.getID(), t);
			}

			if (leaving == 1)
				leavingAtom.add(atom);
		}
	}

	/**
	 * @param bond
	 * @param leaving
	 */
	private void addBondInPseudoMolecule(IBond bond, int leaving) {
		if (!bondsInPseudoMol.contains(bond.getID())) {
			bondsInPseudoMol.add(bond.getID());
			if (leaving == 1) {
				bond = formatBond(bond);
				leavingBond.add(bond);
				pseudoMolecule.addBond(bond);
			}
			else {
				if (repetitionToProcess){
					bond = formatBond(bond);
					pseudoMolecule.addBond(bond);
				}
				else
					pseudoMolecule.addBond(bond);
			}
			//at least one atom has to be in product to be consider as a reaction centre
			if (bond.getProperty(BOND_CHANGE_INFORMATION) != null) {
				if ((int) bond.getProperty(BOND_CHANGE_INFORMATION) == additionalConstants.BOND_FORMED)
					bondFormedList.add(bond);
				if ((int) bond.getProperty(BOND_CHANGE_INFORMATION) == additionalConstants.BOND_CLEAVED)
					bondCleavedList.add(bond);	
				if ((int) bond.getProperty(BOND_CHANGE_INFORMATION) == additionalConstants.BOND_ORDER)
					bondOrderList.add(bond);
				//set first index atom to avoid multiple atom with the same ID in the pseudoMolecule (same reactant interacts multiple times)
				if (bond.getProperty(BOND_CHANGE_INFORMATION) != null) {
					bond.setAtom(indexAtomsPseudoMolUnique.get(bond.getAtom(0).getID()).iterator().next(), 0);
					bond.setAtom(indexAtomsPseudoMolUnique.get(bond.getAtom(1).getID()).iterator().next(), 1);
				}
				if (bond.getProperty(BOND_CHANGE_INFORMATION) != null) {
					bond.getAtom(0).setProperty(additionalConstants.REACTION_CENTER, true);
					bond.getAtom(1).setProperty(additionalConstants.REACTION_CENTER, true);
					reactionCenter.add(bond.getAtom(0));
					reactionCenter.add(bond.getAtom(1));
				}
				
			}
			if (bond.getProperty(additionalConstants.BOND_STEREO) != null) 
				bondStereoList.add(bond);
		}

	}

	/**
	 * @param bond
	 * @return
	 */
	private IBond formatBond(IBond bond) {
		IBond newBond = new Bond();
		newBond.setAtom(indexAtomsPseudoMolUnique.get(bond.getAtom(0).getID()).iterator().next(), 0);
		newBond.setAtom(indexAtomsPseudoMolUnique.get(bond.getAtom(1).getID()).iterator().next(), 1);
		newBond.setElectronCount(bond.getElectronCount());
		newBond.setID(bond.getID());
		newBond.setFlags(bond.getFlags());
		newBond.setIsAromatic(bond.isAromatic());
		newBond.setIsInRing(bond.isInRing());
		newBond.setNotification(bond.getNotification());
		newBond.setOrder(bond.getOrder());
		newBond.setProperties(bond.getProperties());
		newBond.setStereo(bond.getStereo());
		return newBond;
	}

	/**
	 * @return
	 */
	public HashMap<String, Integer> atomRepetition() {
		HashMap<String, Integer> numberOfRepetitions = new HashMap<String, Integer>();
		for (Entry<String, LinkedHashSet<IAtom>> set : indexAtomsPseudoMol.entrySet()) {
			numberOfRepetitions.put(set.getKey(), set.getValue().size());
		}

		return numberOfRepetitions;
	}

	/**
	 * @param mol
	 * @return
	 * @throws CDKException
	 */
	private IAtomContainer cleanMolecule(IAtomContainer mol) throws CDKException {
		//use StructureDiagramGenerator2 because StructureDiagramGenerator modify the stereo and can remove it
		StructureDiagramGenerator2 sdg = new StructureDiagramGenerator2();
		sdg.setMolecule(mol, false);
		sdg.generateCoordinates();

		return sdg.getMolecule();
	}

	/**
	 * @param newPseudoMolecule
	 * @return
	 * @throws CDKException
	 */
	public String GetPseudoMoleculeSmiles(IAtomContainer newPseudoMolecule) throws CDKException {
		SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Default | SmiFlavor.UseAromaticSymbols);
		return sg.create(newPseudoMolecule);
	}

	/**
	 * @param reaction
	 * @throws CDKException
	 */
	public void reactionAnnotator(IReaction reaction) throws CDKException {
		//idofTheAtom:indexOfTheAtom
		Map<Integer,Integer> indexAtomsInReactants = new HashMap<Integer,Integer>();
		Map<Integer,Integer> indexAtomsInProducts = new HashMap<Integer,Integer>();
		//idofTheBond:indexOfTheBond
		Map<String,Integer> indexBondsInReactants = new HashMap<String,Integer>();
		Map<String,Integer> indexBondsInProducts = new HashMap<String,Integer>();
		
		//Map to correct unbalanced reactions of type substitution idOfSubstituteAtom:idsOfAdjacentAtomsInReactant
		Map<Integer,List<Integer>> substitutions = new HashMap<Integer,List<Integer>>();

		//merge all reactant and product to find difference
		IAtomContainer reactants = reaction.getReactants().getAtomContainer(0);
		IAtomContainer products =  reaction.getProducts().getAtomContainer(0);


		for (int i = 1; i < reaction.getProducts().getAtomContainerCount(); i++) {
			IAtomContainer ac = reaction.getProducts().getAtomContainer(i);
			products.add(ac);
		}
		
		//reference atoms index
		for (int j = 0; j < products.getAtomCount(); j++) {
			indexAtomsInProducts.put(products.getAtom(j).getProperty(CDKConstants.ATOM_ATOM_MAPPING), j);
		}
		//reference bonds index
		for (int j = 0; j < products.getBondCount(); j++) {
			IBond bond = products.getBond(j);
			System.out.println(bond.getID());
			indexBondsInProducts.put(bond.getID(), j);
		}

		for (int i = 0; i < reaction.getReactants().getAtomContainerCount(); i++) {
			IAtomContainer reactant = reaction.getReactants().getAtomContainer(i);

			//for (IStereoElement se : ac.stereoElements())
			//	rSE.add(se);

			for (int j = 0; j < reactant.getAtomCount(); j++) {
				IAtom rA = reactant.getAtom(j);
				indexAtomsInReactants.put(rA.getProperty(CDKConstants.ATOM_ATOM_MAPPING), j);
				if (!indexAtomsInProducts.containsKey(rA.getProperty(CDKConstants.ATOM_ATOM_MAPPING))) 
					rA.setProperty(additionalConstants.LEAVING_ATOM, true);
				else {
					IAtom pA = products.getAtom(indexAtomsInProducts.get(rA.getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
					if (rA.getFormalCharge() > pA.getFormalCharge()) {
						pA.setProperty(additionalConstants.CHARGE_CHANGE, true);
						pA.setProperty(additionalConstants.ATOM_CHARGE_CHANGE, additionalConstants.LOOSE_ONE_CHARGE);

					}
					else if (rA.getFormalCharge() < pA.getFormalCharge()) {
						pA.setProperty(additionalConstants.CHARGE_CHANGE, true);
						pA.setProperty(additionalConstants.ATOM_CHARGE_CHANGE, additionalConstants.GAIN_ONE_CHARGE);

					}
				}	
			}
			for (int j = 0; j < reactant.getBondCount(); j++) {
				IBond rB = reactant.getBond(j);
				indexBondsInReactants.put(rB.getID(), j);
				rB.setProperty(additionalConstants.BOND_BOND_MAPPING, j+1);
				if (indexBondsInProducts.containsKey(rB.getID())) {
					IAtom pABegin = products.getAtom(indexAtomsInProducts.get(rB.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
					IAtom pAEnd = products.getAtom(indexAtomsInProducts.get(rB.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
					IBond pB = products.getBond(pABegin,pAEnd);
					pB.setProperty(additionalConstants.BOND_BOND_MAPPING, j+1);
					if (encodeOrder(rB) != encodeOrder(pB) && !(rB.isAromatic() == true && pB.isAromatic() == true)) {
						pB.setProperty(additionalConstants.BOND_ORDER, true);
						pB.setProperty(BOND_CHANGE_INFORMATION, additionalConstants.BOND_ORDER);
						pB.setProperty(additionalConstants.BOND_ORDER_CHANGE, getBondOderChangeProperty(rB.getOrder(), pB.getOrder()));
						pB.getBegin().setProperty(additionalConstants.REACTION_CENTER, true);
						pB.getEnd().setProperty(additionalConstants.REACTION_CENTER, true);
					}
					if (!rB.getStereo().equals(pB.getStereo())) {
						pB.setProperty(additionalConstants.BOND_STEREO, true);
						pB.setProperty(additionalConstants.BOND_STEREO_CHANGE, getBondStereoChangeProperty(rB.getStereo(), pB.getStereo()));
						//pB.getBegin().setProperty(additionalConstants.REACTION_CENTER, true);
						//pB.getEnd().setProperty(additionalConstants.REACTION_CENTER, true);
					}
				}
				else {
					boolean isCleaved = false;
					int beginR = rB.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING);
					int endR = rB.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING);
					List<Integer> adj = new ArrayList<Integer>();
					if ((indexAtomsInProducts.containsKey(beginR) || indexAtomsInProducts.containsKey(endR))) {
						int ref = -1;
						if (indexAtomsInProducts.containsKey(beginR)) {
							ref = beginR;
							if (!indexAtomsInProducts.containsKey(endR)) {
								isCleaved = true;
							}
							else {
								reactant.getConnectedAtomsList(rB.getBegin())
								.forEach(atom -> adj.add(atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
								for (int atom : adj) {
									if (indexAtomsInProducts.containsKey(atom)) {
										isCleaved = true;
									}
								}
							}
						}
						else if (indexAtomsInProducts.containsKey(endR)) {
							ref = endR;
							if (!indexAtomsInProducts.containsKey(beginR)) {
								isCleaved = true;
							}
							else {
								reactant.getConnectedAtomsList(rB.getEnd())
								.forEach(atom -> adj.add(atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
								for (int atom : adj) {
									if (indexAtomsInProducts.containsKey(atom)) {
										isCleaved = true;
									}
								}
							}
							
						}
						if (isCleaved) {
							rB.setProperty(BOND_CHANGE_INFORMATION, additionalConstants.BOND_CLEAVED);
							rB.getBegin().setProperty(additionalConstants.REACTION_CENTER, true);
							rB.getEnd().setProperty(additionalConstants.REACTION_CENTER, true);
						}
						//probable substitution
						else {
							substitutions.put(ref, adj);
						}
					}
					
				}
			}

			List<IStereoElement> rSE = new ArrayList<IStereoElement>();

			for (IStereoElement se : reactant.stereoElements())
				rSE.add(se);

			if (rSE.size() > 0) {
			
				//Annotate  stereoCenter
				for (IStereoElement se : rSE) {
					if (se.getFocus() instanceof IAtom) {
						IAtom atom = (IAtom) se.getFocus();
						atom.setProperty(additionalConstants.IS_STEREOCENTER, true);
						configStereoProperties(se, atom);
					}
					else if (se.getFocus() instanceof IBond) {
						IAtom begin = ((IBond) se.getFocus()).getAtom(0);
						IAtom end = ((IBond) se.getFocus()).getAtom(1);
						begin.setProperty(additionalConstants.IS_STEREOCENTER, true);
						configStereoProperties(se, begin);
						end.setProperty(additionalConstants.IS_STEREOCENTER, true);
						configStereoProperties(se, end);
					}
				}
			}
			
			if (i > 0)
				reactants.add(reactant);
		}

		int countBondMapped = reactants.getBondCount();
		for (int i = 0; i < reaction.getProducts().getAtomContainerCount(); i++) {
			IAtomContainer product = reaction.getProducts().getAtomContainer(i);

			for (int j = 0; j < product.getBondCount(); j++) {
				IBond pB = product.getBond(j);
				if (!indexBondsInReactants.containsKey(pB.getID())) {
//					System.out.println("formed " + pB.getBegin().getSymbol() + pB.getBegin().getID() + " " +
//							pB.getEnd().getSymbol() + pB.getEnd().getID() );
					pB.setProperty(additionalConstants.BOND_BOND_MAPPING, countBondMapped);
					pB.setProperty(BOND_CHANGE_INFORMATION, additionalConstants.BOND_FORMED);
					pB.getBegin().setProperty(additionalConstants.REACTION_CENTER, true);
					pB.getEnd().setProperty(additionalConstants.REACTION_CENTER, true);
					countBondMapped++;
					
					//correct substitution for unbalanced reactions ex CN(C)C=O.COC1=CC(=C(C=C1)C(O)=O)[N+]([O-])=O>>COC1=CC(=C(C=C1)C(N)=O)[N+]([O-])=O
//					if (!substitutions.isEmpty()) {
//						int beginP = pB.getBegin().getProperty(CDKConstants.ATOM_ATOM_MAPPING);
//						int endP = pB.getEnd().getProperty(CDKConstants.ATOM_ATOM_MAPPING);
//						
//						if (substitutions.containsKey(beginP)) {
//							IAtom commonAtomInReactantAndProduct = reactants.getAtom(indexAtomsInReactants.get(endP));
//							//get adj atoms of other atom (which is in reactant too;
//							List<IBond> con = reactants.getConnectedBondsList(commonAtomInReactantAndProduct);
//							for (IBond bond : con) {
//								IAtom uniqueAtomInReactant = bond.getOther(commonAtomInReactantAndProduct);
//								if (uniqueAtomInReactant.getProperty(additionalConstants.REACTION_CENTER) != null) {
//									if ((int)bond.getProperty(BOND_CHANGE_INFORMATION) == additionalConstants.BOND_CLEAVED) {
//										//make formed bond with adj
//									}
//								}
//							}
//						}
//					}
				}
			}

			List<IStereoElement> pSE = new ArrayList<IStereoElement>();

			for (IStereoElement se : product.stereoElements())
				pSE.add(se);

			//Annotate  stereoCenter
			if (pSE.size() > 0) {
			
				//Annotate  stereoCenter
				for (IStereoElement se : pSE) {
					if (se.getFocus() instanceof IAtom) {
						IAtom atom = (IAtom) se.getFocus();
						atom.setProperty(additionalConstants.IS_STEREOCENTER, true);
						configStereoProperties(se, atom);
					}
					else if (se.getFocus() instanceof IBond) {
						IAtom begin = ((IBond) se.getFocus()).getAtom(0);
						IAtom end = ((IBond) se.getFocus()).getAtom(1);
						begin.setProperty(additionalConstants.IS_STEREOCENTER, true);
						configStereoProperties(se, begin);
						end.setProperty(additionalConstants.IS_STEREOCENTER, true);
						configStereoProperties(se, end);
					}
				}
			}
		}
	}
	
	
	/**
	 * LEFT == OPPOSITE == ANTI_CLOCKWISE == 1
	 * RIGHT == TOGETHER == CLOCKWISE == 2
	 * @param se
	 * @param a
	 */
	private void configStereoProperties(IStereoElement se, IAtom a) {
		String symbol = null;
		String type = null;
		String shorthand = null;
		if (se instanceof ITetrahedralChirality || se instanceof TetrahedralChirality) {
			ITetrahedralChirality tc = (ITetrahedralChirality) se;
			ITetrahedralChirality.Stereo stereo = tc.getStereo();
			if (stereo.equals(ITetrahedralChirality.Stereo.CLOCKWISE)) {
				symbol = "@TH2";
				type = "Tetrahedral";
				shorthand = "CLOCKWISE";
			}
			else if (stereo.equals(ITetrahedralChirality.Stereo.ANTI_CLOCKWISE)) {
				symbol = "@TH1";
				type = "Tetrahedral";
				shorthand = "ANTI_CLOCKWISE";
			}
		} 
		else if (se instanceof ExtendedTetrahedral) {
			ExtendedTetrahedral et = (ExtendedTetrahedral) se;
			ITetrahedralChirality.Stereo stereo = et.winding();
			if (stereo.equals(ITetrahedralChirality.Stereo.CLOCKWISE)) {
				symbol = "@AL2";
				type = "Tetrahedral";
				shorthand = "CLOCKWISE";
			}
			else if (stereo.equals(ITetrahedralChirality.Stereo.ANTI_CLOCKWISE)) {
				symbol = "@AL1";
				type = "Tetrahedral";
				shorthand = "ANTI_CLOCKWISE";
			}
		} 
		else if (se instanceof SquarePlanar) {
			SquarePlanar sp = (SquarePlanar) se;
			int configOrder = sp.getConfigOrder();
			if (configOrder == 1) {
				symbol = "@SP1";
				type = "SquarePlanar";
				shorthand = "";
			}
			else if (configOrder == 2) {
				symbol = "@SP2";
				type = "SquarePlanar";
				shorthand = "";
			}
			else if (configOrder == 3) {
				symbol = "@SP3";
				type = "SquarePlanar";
				shorthand = "";
			}
		}
		else if (se instanceof TrigonalBipyramidal) {
			TrigonalBipyramidal tb = (TrigonalBipyramidal) se;
			int configOrder = tb.getConfigOrder();
			if (configOrder == IStereoElement.RIGHT) {
				symbol = "@TB2";
				type = "TrigonalBipyramidal";
				shorthand = "CLOCKWISE";
			}
			else if (configOrder == IStereoElement.LEFT) {
				symbol = "@TB1";
				type = "TrigonalBipyramidal";
				shorthand = "ANTI_CLOCKWISE";
			}
		}
		else if (se instanceof Octahedral) {
			Octahedral oh = (Octahedral) se;
			int configOrder = oh.getConfigOrder();
			if (configOrder == IStereoElement.RIGHT) {
				symbol = "@OH2";
				type = "Octahedral";
				shorthand = "CLOCKWISE";
			}
			else if (configOrder == IStereoElement.LEFT) {
				symbol = "@OH1";
				type = "Octahedral";
				shorthand = "ANTI_CLOCKWISE";
			}
		}
		else if (se instanceof IDoubleBondStereochemistry || se instanceof DoubleBondStereochemistry) {
			IDoubleBondStereochemistry db = (IDoubleBondStereochemistry) se;
			int config = db.getConfigOrder();
			if (config == IStereoElement.TOGETHER) {
				symbol = "@DB2";
				type = "DoubleBondStereochemistry";
				shorthand = "CLOCKWISE";
			}
			else if (config == IStereoElement.OPPOSITE) {
				symbol = "@DB1";
				type = "DoubleBondStereochemistry";
				shorthand = "ANTI_CLOCKWISE";
			}
		}
		else if (se instanceof ExtendedCisTrans) {
			ExtendedCisTrans ct = (ExtendedCisTrans) se;
			int config = ct.getConfigOrder();
			if (config == IStereoElement.TOGETHER) {
				symbol = "@CT2";
				type = "ExtendedCisTrans";
				shorthand = "CLOCKWISE";
			}
			else if (config == IStereoElement.OPPOSITE) {
				symbol = "@CT1";
				type = "ExtendedCisTrans";
				shorthand = "ANTI_CLOCKWISE";

			}
			else if (se instanceof Atropisomeric) {
				Atropisomeric ap = (Atropisomeric) se;
				int configOrder = ap.getConfigOrder();
				if (configOrder == IStereoElement.RIGHT) {
					symbol = "@AP2";
					type = "Atropisomeric";
					shorthand = "CLOCKWISE";
				}
				else if (configOrder == IStereoElement.LEFT) {
					symbol = "@AP1";
					type = "Atropisomeric";
					shorthand = "ANTI_CLOCKWISE";
				}
			}
		}
		a.setProperty(additionalConstants.STEREO_TYPE, type);
		a.setProperty(additionalConstants.STEREO_SYMBOL, symbol);
		a.setProperty(additionalConstants.STEREO_SHORTHAND, shorthand);
	}
	
	
	/**
	 * @param bond
	 * @return
	 */
	private int encodeOrder(IBond bond) {
		if (bond.isAromatic()) {
			return 9;
		}
		else {
			return bond.getOrder().numeric();
		}
	}
	
	/**
	 * Return bond order change between the reactant and the product
	 * @param r
	 * @param p
	 * @return
	 */
	private int getBondOderChangeProperty(IBond.Order r, IBond.Order p) {
		if (r.numeric() > p.numeric()) 
			return additionalConstants.BOND_ORDER_REDUCED;
		else if (r.numeric() < p.numeric()) 
			return additionalConstants.BOND_ORDER_GAIN;
		else
			return -1;
	}
	

	/**
	 * @param r
	 * @param p
	 * @return
	 */
	private int getBondStereoChangeProperty(IBond.Stereo r, IBond.Stereo p) {
		if (r == null && p != null) {
			if (p.equals(IBond.Stereo.UP)) 
				return additionalConstants.NONE_TO_UP;
			else if (p.equals(IBond.Stereo.DOWN)) 
				return additionalConstants.NONE_TO_DOWN;
			else if (p.equals(IBond.Stereo.UP_INVERTED)) 
				return additionalConstants.NONE_TO_UP_INVERTED;
			else if (p.equals(IBond.Stereo.DOWN_INVERTED)) 
				return additionalConstants.NONE_TO_DOWN_INVERTED;
			else if (p.equals(IBond.Stereo.Z)) 
				return additionalConstants.NONE_TO_Z;
			else if (p.equals(IBond.Stereo.E)) 
				return additionalConstants.NONE_TO_E;
			else
				return -1;
		}
		else if (r != null && p == null)  {
			if (r.equals(IBond.Stereo.NONE)) 
				return additionalConstants.UP_TO_NONE;
			else if (r.equals(IBond.Stereo.NONE)) 
				return additionalConstants.DOWN_TO_NONE;
			else if (r.equals(IBond.Stereo.NONE)) 
				return additionalConstants.UP_INVERTED_TO_NONE;
			else if (r.equals(IBond.Stereo.NONE)) 
				return additionalConstants.DOWN_INVERTED_TO_NONE;
			else if (r.equals(IBond.Stereo.NONE)) 
				return additionalConstants.Z_TO_NONE;
			else if (r.equals(IBond.Stereo.NONE)) 
				return additionalConstants.E_TO_NONE;
			else
				return -1;
		}
		else {
			if (r.equals(IBond.Stereo.DOWN) && p.equals(IBond.Stereo.UP)) 
				return additionalConstants.DOWN_TO_UP;
			else if (r.equals(IBond.Stereo.UP) && p.equals(IBond.Stereo.DOWN)) 
				return additionalConstants.UP_TO_DOWN;
			else if (r.equals(IBond.Stereo.DOWN_INVERTED) && p.equals(IBond.Stereo.UP_INVERTED)) 
				return additionalConstants.DOWN_INVERTED_TO_UP_INVERTED;
			else if (r.equals(IBond.Stereo.UP_INVERTED) && p.equals(IBond.Stereo.DOWN_INVERTED)) 
				return additionalConstants.UP_INVERTED_TO_DOWN_INVERTED;
			else if (r.equals(IBond.Stereo.E) && p.equals(IBond.Stereo.Z)) 
				return additionalConstants.E_TO_Z;
			else if (r.equals(IBond.Stereo.Z) && p.equals(IBond.Stereo.E)) 
				return additionalConstants.Z_TO_E;
			else if (r.equals(IBond.Stereo.NONE) && p.equals(IBond.Stereo.UP)) 
				return additionalConstants.NONE_TO_UP;
			else if (r.equals(IBond.Stereo.NONE) && p.equals(IBond.Stereo.DOWN)) 
				return additionalConstants.NONE_TO_DOWN;
			else if (r.equals(IBond.Stereo.NONE) && p.equals(IBond.Stereo.UP_INVERTED)) 
				return additionalConstants.NONE_TO_UP_INVERTED;
			else if (r.equals(IBond.Stereo.NONE) && p.equals(IBond.Stereo.DOWN_INVERTED)) 
				return additionalConstants.NONE_TO_DOWN_INVERTED;
			else if (r.equals(IBond.Stereo.NONE) && p.equals(IBond.Stereo.Z)) 
				return additionalConstants.NONE_TO_Z;
			else if (r.equals(IBond.Stereo.NONE) && p.equals(IBond.Stereo.E)) 
				return additionalConstants.NONE_TO_E;
			else if (r.equals(IBond.Stereo.DOWN) && p.equals(IBond.Stereo.UP_INVERTED)) 
				return additionalConstants.DOWN__TO_UP_INVERTED;
			else if (r.equals(IBond.Stereo.UP) && p.equals(IBond.Stereo.DOWN_INVERTED)) 
				return additionalConstants.UP_TO_DOWN_INVERTED;
			else if (r.equals(IBond.Stereo.DOWN_INVERTED) && p.equals(IBond.Stereo.UP)) 
				return additionalConstants.DOWN_INVERTED_TO_UP;
			else if (r.equals(IBond.Stereo.UP_INVERTED) && p.equals(IBond.Stereo.DOWN)) 
				return additionalConstants.UP_INVERTED_TO_DOWN;
			else
				return -1;
		}
	}
	
	

	/**
	 * Use atom in products to build an index. If an atom is repeated (multiple atoms with the same atom atom mapping)
	 * , repetitionToProcess is set true and duplicate atom will be deleted to avoid duplicate in the pseudoMolecules
	 * ex: 2A + 1B -> 2AB  The molecule A is to time in the product and has to be considered as the same molecule repeated 2 times
	 * (cf SI figure Stoichiometry management in publication)
	 * @param atom
	 * @return
	 */
	private IAtomContainer buildIndexInPseudoMolecule(IAtomContainerSet set) {
		IAtomContainer aggregate = DefaultChemObjectBuilder.getInstance().newAtomContainer();
		for (IAtomContainer ac : set.atomContainers()) {
			aggregate.add(ac);
			for (IAtom atom : ac.atoms()) {
				// index Atom to identify if an atom is repeated in product 
				if (indexAtomsPseudoMolUnique.get(atom.getID()) == null) {
					repetitionToProcess = true;
					LinkedHashSet<IAtom> t = new LinkedHashSet<IAtom>();
					t.add(atom);
					indexAtomsPseudoMolUnique.put(atom.getID(), t);
				} else {
					LinkedHashSet<IAtom> t = indexAtomsPseudoMolUnique.get(atom.getID());
					t.add(atom);
					indexAtomsPseudoMolUnique.put(atom.getID(), t);
				}
			}
		}
		indexAtomsPseudoMol = new HashMap<String,LinkedHashSet<IAtom>>(indexAtomsPseudoMolUnique);
		return aggregate;
	}
	

	/**
	 * Remove smartly duplicate atoms. A BFS is looking for neighbors in order to keep a whole fragment
	 * @param ac
	 */
	private void duplicateAtomManagement(IAtomContainer ac) {
		for (Entry<String,LinkedHashSet<IAtom>> e : indexAtomsPseudoMol.entrySet()) {
			String id = e.getKey();
			LinkedHashSet<IAtom> atoms = indexAtomsPseudoMolUnique.get(id);
			if (atoms.size() > 1) {
				Set<IAtom> sphere = new HashSet<IAtom>();
				sphere.add(atoms.iterator().next());
				List<IAtom> selected  = getAllConnectedDuplicateAtoms(sphere, ac, 
						new ArrayList<IAtom>(), new HashSet<IAtom>());
				selected.add(atoms.iterator().next());
				//define the atoms in selected as the reference atoms to be kept
				for (IAtom atom : selected) {
					LinkedHashSet<IAtom> temp = new LinkedHashSet<IAtom>();
					temp.add(atom);
					indexAtomsPseudoMolUnique.put(atom.getID(), temp);
				}
			}
		}
	}
	
	/**
	 * @param sphere
	 * @param atomContainer
	 * @param atomList
	 * @param visited
	 * @return
	 */
	private List<IAtom> getAllConnectedDuplicateAtoms(Set<IAtom> sphere, IAtomContainer atomContainer, 
			List<IAtom> atomList, Set<IAtom> visited) {
		IAtom nextAtom;
		Set<IAtom> newSphere = new HashSet<IAtom>();

		for (IAtom atom : sphere) {
			List<IBond> connectedBonds = atomContainer.getConnectedBondsList(atom);
			for (IBond bond : connectedBonds) {
				nextAtom = bond.getOther(atom);
				//only add non visited and duplicate atoms
				if (!visited.contains(nextAtom)) {
					if (!sphere.contains(nextAtom)  && indexAtomsPseudoMolUnique.get(nextAtom.getID()).size() > 1) {
						newSphere.add(nextAtom);
						if (!atomList.contains(nextAtom)) atomList.add(nextAtom);
					}
					visited.add(nextAtom);
				}
			}
		}
		if (newSphere.size() > 0) {
			getAllConnectedDuplicateAtoms(newSphere, atomContainer, atomList, visited);
		}
		return atomList;
	}
	
	/**
	 * @param name
	 * @param newPseudoMolecule
	 * @throws IOException
	 * @throws CDKException
	 * @throws org.openscience.cdk.exception.CDKException 
	 */
	public void writePseudoMoleculeSDFile(String name, String path, IAtomContainer newPseudoMolecule) throws IOException, CDKException, org.openscience.cdk.exception.CDKException {
		File f = new File(new File(path).getCanonicalPath() + separator + name + ".sd");
		try (MDLV2000Writer writer = new MDLV2000Writer(new FileWriter(f))) {
			writer.write(newPseudoMolecule);
			writer.close();
		}
	}

	/**
	 * @param parentPath
	 * @param name
	 * @throws IOException
	 * @throws CDKException
	 * @throws CloneNotSupportedException
	 */
	public void writePseudoMoleculeImage(String parentPath, String name) throws IOException, CDKException, CloneNotSupportedException {
		new org.openscience.cdk.depict.DepictionGenerator().withHighlight(getBondOrderList(), Color.GREEN)
		.withHighlight(getBondFormedList(), Color.BLUE)
		.withHighlight(getBondCleavedList(), Color.RED)
		.withHighlight(getBondOrderList(), Color.GREEN)
		.withHighlight(getReactioncenter(), Color.LIGHT_GRAY)
		.withOuterGlowHighlight().withAtomColors().withAtomMapNumbers()
		.depict(new AtomContainer(pseudoMolecule)).writeTo(new File(parentPath, name+"_pseudoMolecule.pdf").toString());
	}
	
	/**
	 * @return
	 */
	public Set<IBond> getBondOrderList() {
		return bondOrderList;
	}

	/**
	 * @return
	 */
	public Set<IBond> getBondStereoList() {
		return bondStereoList;
	}

	/**
	 * @return
	 */
	public Set<IBond> getBondFormedList() {
		return bondFormedList;
	}

	/**
	 * @return
	 */
	public Set<IBond> getBondCleavedList() {
		return bondCleavedList;
	}

	/**
	 * @return
	 */
	public Set<IAtom> getReactioncenter() {
		return reactionCenter;
	}

	/**
	 * @return
	 */
	public Map<String, LinkedHashSet<IAtom>> getIndexAtomsPseudoMolecule() {
		return indexAtomsPseudoMolUnique;
	}

	/**
	 * @param indexAtomsPseudoMolecule
	 */
	public void setIndexAtomsPseudoMolecule(HashMap<String, LinkedHashSet<IAtom>> indexAtomsPseudoMolecule) {
		this.indexAtomsPseudoMolUnique = indexAtomsPseudoMolecule;
	}


}
