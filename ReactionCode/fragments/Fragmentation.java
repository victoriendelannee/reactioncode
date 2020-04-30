package com.nih.fragments;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IStereoElement;

import com.nih.util.tools;

public class Fragmentation {

	private final String IS_SIDECHAIN_ATOM    = "sidechainatom";
	private final String IS_LINKER_ATOM       = "linkeratom";
	private final String IS_CONNECTED_TO_RING = "rcon";
	private final String IS_IN_MURCKO_SCAFFOLD    = "murcko";
	private final String IS_SIDECHAIN    = "sidechain";
	private final String IS_LINKER       = "linker";
	
	private boolean HAS_LINKER = false;

	List<IAtomContainer> linkers;
	List<IAtomContainer> sideChains;
	List<IAtomContainer> ringSystems;
	
	//Connections between a ring system (key) and their linker (value)
	Map<IAtomContainer,List<IAtomContainer>> ringSystemLinkerConnections;
	//Connections between a ring system (key) and their side chains and linker (value)
	Map<IAtomContainer,List<IAtomContainer>> ringSystemSideChainConnections;
	//Atom Connection between ringSystem (key) and a linker/side chain (value)
	Map<String,String> connectedAtoms;
	//Connections between a ring system (key) and their linkers and side chains (value)
	Map<IAtomContainer,List<IAtomContainer>> connections = new HashMap<IAtomContainer,List<IAtomContainer>>();
	
	//key = ring atom ID : value = bond with the ring atom and the connector atom (*), which will have to be connected to an aliphatic chain
	Map<String,IBond> ringToConnector;




	//TODO Find a way to Manage stereo 
	/**
	 * Split the side chains and the linkers from the ring systems and set IAtomContainerSet containing the <\br>
	 * rings, the side chains and the linker. The ring systems are connected to the sides chains and linkers by <\br>
	 * bonds, where atoms a1 ISINRING and a2 is not. For all ring system, these bonds are kept but all atoms a2 <\br>
	 * are replaced by a generic group (L for a linker and S for a side chain) following by a number (to later reference it in the database). For all <\br>
	 * linkers and side chains, the atom a1 is replace by a generic R (for ring) atom following by number. this atom has the <\br>
	 * number than the one in the bond to connect in the ring (ex: S1 <-> R1).
	 * @param ac
	 * @return
	 * @throws CloneNotSupportedException 
	 */
	public void splitLinkersSideChainsAndRingSystem(IAtomContainer ac) throws CloneNotSupportedException {
		linkers = new ArrayList<IAtomContainer>();
		sideChains = new ArrayList<IAtomContainer>();
		ringSystems = new ArrayList<IAtomContainer>();
		
		ringSystemLinkerConnections = new HashMap<IAtomContainer,List<IAtomContainer>>();
		ringSystemSideChainConnections = new HashMap<IAtomContainer,List<IAtomContainer>>();
		connectedAtoms = new HashMap<String,String>();
		ringToConnector = new HashMap<String,IBond>();
		
		Set<IBond> bondToRemove = new HashSet<IBond>();
		Set<IBond> bondToAdd = new HashSet<IBond>();
		int counter = 1;

		for (IBond b : ac.bonds()) {
			IBond bRing;
			IBond bChain;
			Atom other1;
			Atom other2;
			if (!b.getAtom(0).getFlag(CDKConstants.ISINRING) && b.getAtom(1).getFlag(CDKConstants.ISINRING)) {
				if ((boolean) b.getAtom(0).getProperty(IS_LINKER_ATOM)) {
					//IS_CONNECTED_TO_RING : ISINRING
//					Map<IAtom,IAtom> linkerAtoms = new HashMap<IAtom,IAtom>();
//
//					Set<IAtom> sphere = new HashSet<IAtom>();
//					sphere.add(b.getAtom(0));
//					getOtherAtomsConnectedToARingSystem(sphere, ac, linkerAtoms);

//					for (Entry<IAtom,IAtom> e : linkerAtoms.entrySet()) {
//						IAtom connectedToRing = e.getKey();
//						IAtom inRing = e.getValue();
						IAtom connectedToRing = b.getAtom(0);
						IAtom inRing = b.getAtom(1);
						IBond link = ac.getBond(connectedToRing, inRing);
						if (link != null) {
							bRing = new Bond();
							other1 = new Atom();
							other1.setSymbol("*");
							other1.setProperty("connection", "L_"+counter);
							other1.setID("L_"+counter);
							bRing.setAtom(inRing, 0);
							bRing.setAtom(other1, 1);
							bRing.setOrder(b.getOrder());
							bRing.setProperty("connector", true);
//							bRing.setProperty("linkerConnector", true);
							System.out.println(bRing);

							bChain = new Bond();
							other2 = new Atom();
							other2.setSymbol("*");
							other2.setProperty("connection", "R_"+counter);
							other2.setID("R_"+counter);
							other2.setProperty("isLinkerConnector", true);
							bChain.setAtom(connectedToRing, 0);
							bChain.setAtom(other2, 1);
							bChain.setOrder(b.getOrder());
							bChain.setProperty("connector", true);
							bChain.setProperty("isLinkerConnector", true);
							
							ac.addAtom(other1);
							ac.addAtom(other2);
							
							bondToAdd.add(bChain);
							bondToAdd.add(bRing);
							
							ringToConnector.put(inRing.getID(), bRing);
							connectedAtoms.put(other1.getID(), other2.getID());
							
							bondToRemove.add(b);
							counter++;
//						}

					}
				}
				else {
					bRing = new Bond();
					other1 = new Atom();
					other1.setSymbol("*");
					other1.setProperty("connection", "S_"+counter);
					other1.setID("S_"+counter);
					bRing.setAtom(b.getAtom(1), 0);
					bRing.setAtom(other1, 1);
					bRing.setOrder(b.getOrder());
					bRing.setProperty("connector", true);
//					bRing.setProperty("sideConnector", true);

					bChain = new Bond();
					other2 = new Atom();
					other2.setSymbol("*");
					other2.setProperty("connection", "R_"+counter);
					other2.setID("R_"+counter);
					other2.setProperty("isSideConnector", true);
					bChain.setAtom(b.getAtom(0), 0);
					bChain.setAtom(other2, 1);
					bChain.setOrder(b.getOrder());
					bChain.setProperty("connector", true);
					bChain.setProperty("isSideConnector", true);
					
					ac.addAtom(other1);
					ac.addAtom(other2);
					
					bondToAdd.add(bChain);
					bondToAdd.add(bRing);
					
					ringToConnector.put(b.getAtom(1).getID(), bRing);
					connectedAtoms.put(other1.getID(), other2.getID());
					
					bondToRemove.add(b);
					counter++;
				}
			}
			else if (!b.getAtom(1).getFlag(CDKConstants.ISINRING) && b.getAtom(0).getFlag(CDKConstants.ISINRING)) {
				if ((boolean) b.getAtom(1).getProperty(IS_LINKER_ATOM)) {
					//IS_CONNECTED_TO_RING : ISINRING
//					Map<IAtom,IAtom> linkerAtoms = new HashMap<IAtom,IAtom>();
//
//					Set<IAtom> sphere = new HashSet<IAtom>();
//					sphere.add(b.getAtom(1));
//					getOtherAtomsConnectedToARingSystem(sphere, ac, linkerAtoms);
//
//					for (Entry<IAtom,IAtom> e : linkerAtoms.entrySet()) {
//						IAtom connectedToRing = e.getKey();
//						IAtom inRing = e.getValue();
					IAtom connectedToRing = b.getAtom(1);
					IAtom inRing = b.getAtom(0);
						IBond link = ac.getBond(connectedToRing, inRing);
						if (link != null) {
							bRing = new Bond();
							other1 = new Atom();
							other1.setSymbol("*");
							other1.setProperty("connection", "L_"+counter);
							other1.setID("L_"+counter);
							bRing.setAtom(inRing, 0);
							bRing.setAtom(other1, 1);
							bRing.setOrder(b.getOrder());
							bRing.setProperty("connector", true);
//							bRing.setProperty("linkerConnector", true);
							System.out.println(bRing);

							bChain = new Bond();
							other2 = new Atom();
							other2.setSymbol("*");
							other2.setProperty("connection", "R_"+counter);
							other2.setID("R_"+counter);
							other2.setProperty("isLinkerConnector", true);
							bChain.setAtom(connectedToRing, 0);
							bChain.setAtom(other2, 1);
							bChain.setOrder(b.getOrder());
							bChain.setProperty("connector", true);
							bChain.setProperty("isLinkerConnector", true);
							
							ac.addAtom(other1);
							ac.addAtom(other2);
							
							bondToAdd.add(bChain);
							bondToAdd.add(bRing);
							
							ringToConnector.put(inRing.getID(), bRing);
							connectedAtoms.put(other1.getID(), other2.getID());
							
							bondToRemove.add(b);
							counter++;
//						}

					}
				}
				else {
					bRing = new Bond();
					other1 = new Atom();
					other1.setSymbol("*");
					other1.setProperty("connection", "S_"+counter);
					other1.setID("S_"+counter);
					bRing.setAtom(b.getAtom(0), 0);
					bRing.setAtom(other1, 1);
					bRing.setOrder(b.getOrder());
					bRing.setProperty("connector", true);
//					bRing.setProperty("sideConnector", true);

					bChain = new Bond();
					other2 = new Atom();
					other2.setSymbol("*");
					other2.setProperty("connection", "R_"+counter);
					other2.setID("R_"+counter);
					other2.setProperty("isSideConnector", true);
					bChain.setAtom(b.getAtom(1), 0);
					bChain.setAtom(other2, 1);
					bChain.setOrder(b.getOrder());
					bChain.setProperty("connector", true);
					bChain.setProperty("isSideConnector", true);
					
					ac.addAtom(other1);
					ac.addAtom(other2);
					
					bondToAdd.add(bChain);
					bondToAdd.add(bRing);
					
					ringToConnector.put(b.getAtom(0).getID(), bRing);
					connectedAtoms.put(other1.getID(), other2.getID());
					
					bondToRemove.add(b);
					counter++;
				}
			}
			else if (!b.getFlag(CDKConstants.ISINRING) && b.getAtom(1).getFlag(CDKConstants.ISINRING) 
					&& b.getAtom(0).getFlag(CDKConstants.ISINRING)) {
				if ((boolean) b.getAtom(0).getProperty(IS_LINKER_ATOM) && (boolean) b.getAtom(1).getProperty(IS_LINKER_ATOM)) {
					IAtom inRing1 = b.getAtom(0);
					IAtom inRing2 = b.getAtom(1);
					IBond bRing1 = new Bond();
					other1 = new Atom();
					other1.setSymbol("*");
					other1.setProperty("connection", "L_"+counter);
					other1.setProperty("connection2", "L0_"+counter);
					other1.setProperty("isLinkerConnector", true);
					other1.setID("L_"+counter);
					bRing1.setAtom(inRing1, 0);
					bRing1.setAtom(other1, 1);
					bRing1.setOrder(b.getOrder());
					bRing1.setProperty("connector", true);
					//bRing.setProperty("linkerConnector", true);

					IBond bRing2 = new Bond();
					other2 = new Atom();
					other2.setSymbol("*");
					other2.setProperty("connection", "L_"+counter);
					other2.setProperty("connection2", "L1_"+counter);
					other2.setProperty("isLinkerConnector", true);
					other2.setID("L_"+counter);
					bRing2.setAtom(inRing2, 0);
					bRing2.setAtom(other2, 1);
					bRing2.setOrder(b.getOrder());
					bRing2.setProperty("connector", true);
					
					IBond linker = new Bond();
					linker.setAtom(other1, 0);
					linker.setAtom(other2, 1);
					linker.setOrder(b.getOrder());
					
					IAtomContainer linkerAc = DefaultChemObjectBuilder.getInstance().newAtomContainer();
					linkerAc.addAtom(other1);
					linkerAc.addAtom(other2);
					linkerAc.addBond(linker);
					if (!linkers.contains(linkerAc))
						linkers.add(linkerAc);

					ac.addAtom(other1);
					ac.addAtom(other2);

					bondToAdd.add(bRing1);
					bondToAdd.add(bRing2);

					ringToConnector.put(inRing1.getID(), bRing1);
					ringToConnector.put(inRing2.getID(), bRing2);
					connectedAtoms.put(other1.getID(), other2.getID());
					connectedAtoms.put(other2.getID(), other1.getID());

					bondToRemove.add(b);
					counter++;

				}
			}
		}
		//add bonds
		for (IBond b : bondToAdd)
			ac.addBond(b);
		//remove Bonds
		for (IBond b : bondToRemove)
			ac.removeBond(b);
		
		ringSystemLinkerConnections = new HashMap<IAtomContainer,List<IAtomContainer>>();
		ringSystemSideChainConnections = new HashMap<IAtomContainer,List<IAtomContainer>>();

		//put IAtomContainer in their right lists: linkers, sideChains, ringSystems;
		IAtomContainerSet set = FragmentUtils.makeAtomContainerSet(ac);
		
		//get index container containing a connection type
		Map<String,Integer> index = new HashMap<String,Integer>();
		
		for (int i = 0; i < set.getAtomContainerCount(); i++) {
			for (IAtom a : set.getAtomContainer(i).atoms()) {
				if (a.getProperty("connection") != null)
					index.put(a.getProperty("connection"), i);
				if (a.getProperty("connection2") != null)
					index.put(a.getProperty("connection2"), i);
			}
		}

		for (IAtomContainer ac2 : set.atomContainers()) {
			boolean checked = false;
			for (IAtom a : ac2.atoms()) {
				if (a.getFlag(CDKConstants.ISINRING) && checked == false){
					ringSystems.add(ac2);
					checked = true;
				}
				if (a.getProperty("isLinkerConnector") != null) {
					System.out.println("ac "+ac2);
					if (a.getProperty("connection").toString().contains("R_")) {
						IAtomContainer query = set.getAtomContainer(index.get("L_"+a.getProperty("connection").toString().replace("R_","")));
						//case where a non ring bond link 2 ring system
						if (query == null) {
							if (a.getProperty("connection2").toString().contains("L0_"))
								query = set.getAtomContainer(index.get("L1_"+a.getProperty("connection").toString().replace("L0_","")));
							else if (a.getProperty("connection2").toString().contains("L1_"))
								query = set.getAtomContainer(index.get("L0_"+a.getProperty("connection").toString().replace("L1_","")));
						}
						ac2.setProperty(IS_LINKER, true);
						if (ringSystemLinkerConnections.containsKey(query)) {
							List<IAtomContainer> list = ringSystemLinkerConnections.get(query);
							list.add(ac2);
							ringSystemLinkerConnections.put(query, list);
						}
						else {
							List<IAtomContainer> list = new ArrayList<IAtomContainer>();
							list.add(ac2);
							ringSystemLinkerConnections.put(query, list);
						}
					}
					if (!linkers.contains(ac2))
						linkers.add(ac2);
					continue;
				}
				else if (a.getProperty("isSideConnector") != null) {
					if (a.getProperty("connection").toString().contains("R_")) {
						IAtomContainer query = set.getAtomContainer(index.get("S_"+a.getProperty("connection").toString().replace("R_","")));
						ac2.setProperty(IS_SIDECHAIN, true);
						if (ringSystemSideChainConnections.containsKey(query)) {
							List<IAtomContainer> list = ringSystemSideChainConnections.get(query);
							list.add(ac2);
							ringSystemSideChainConnections.put(query, list);
						}
						else {
							List<IAtomContainer> list = new ArrayList<IAtomContainer>();
							list.add(ac2);
							ringSystemSideChainConnections.put(query, list);
						}
						sideChains.add(ac2);
						continue;
					}
				}
			}
		}
	}

	//TODO Find a way to Manage stereo 
	/**
	 * Split the side chains and the linkers from the ring systems and set the List containing all ring system <\br>
	 * and make map containing a ringSystem as key and the list of its side chains and linkers connected to it
	 * @param ac
	 * @return
	 * @throws CloneNotSupportedException 
	 */
	public void splitLinkersSideChainsAndRingSystem2(IAtomContainer ac) throws CloneNotSupportedException {
		ringSystems = new ArrayList<IAtomContainer>();
		connections = new HashMap<IAtomContainer,List<IAtomContainer>>();

		Set<IBond> bondToRemove = new HashSet<IBond>();
		Set<IBond> bondToAdd = new HashSet<IBond>();
		int counter = 1;

		for (IBond b : ac.bonds()) {
			IBond bChain;
			Atom other1;
			Atom other2;
			if (!b.getFlag(CDKConstants.ISINRING) && !b.getAtom(0).getFlag(CDKConstants.ISINRING) && b.getAtom(1).getFlag(CDKConstants.ISINRING)) {
				IAtom connectedToRing = b.getAtom(0);
				IAtom inRing = b.getAtom(1);
				inRing.setProperty("connectorID", "R_"+counter);
				IBond link = ac.getBond(connectedToRing, inRing);
				if (link != null) {

					bChain = new Bond();
					other1 = new Atom();
					other1.setSymbol("*");
					other1.setProperty("connection", inRing);
					other1.setID("C_"+counter);
					other1.setProperty("connectorID", "C_"+counter);
					if ((boolean) connectedToRing.getProperty(IS_LINKER_ATOM) == true) {
						other1.setProperty(IS_SIDECHAIN_ATOM, false);
						other1.setProperty(IS_LINKER_ATOM, true);
						HAS_LINKER = true;
					}
					else {
						other1.setProperty(IS_LINKER_ATOM, false);
						other1.setProperty(IS_SIDECHAIN_ATOM, true);
					}
					bChain.setAtom(connectedToRing, 0);
					bChain.setAtom(other1, 1);
					bChain.setOrder(b.getOrder());
//					if ((boolean) connectedToRing.getProperty(IS_LINKER_ATOM) == true)
//						bChain.setProperty(IS_LINKER, true);
//					else
//						bChain.setProperty(IS_SIDECHAIN, true);

					ac.addAtom(other1);

					bondToAdd.add(bChain);

					bondToRemove.add(b);
					counter++;
				}
			}
			else if (!b.getFlag(CDKConstants.ISINRING) && !b.getAtom(1).getFlag(CDKConstants.ISINRING) && b.getAtom(0).getFlag(CDKConstants.ISINRING)) {
				IAtom connectedToRing = b.getAtom(1);
				IAtom inRing = b.getAtom(0);
				inRing.setProperty("connectorID", "R_"+counter);
				IBond link = ac.getBond(connectedToRing, inRing);
				if (link != null) {
					bChain = new Bond();
					other1 = new Atom();
					other1.setSymbol("*");
					other1.setProperty("connection", inRing);
					other1.setID("C_"+counter);
					other1.setProperty("connectorID", "C_"+counter);
					if ((boolean) connectedToRing.getProperty(IS_LINKER_ATOM) == true) {
						other1.setProperty(IS_SIDECHAIN_ATOM, false);
						other1.setProperty(IS_LINKER_ATOM, true);
						HAS_LINKER = true;
					}
					else {
						other1.setProperty(IS_SIDECHAIN_ATOM, true);
						other1.setProperty(IS_LINKER_ATOM, false);
					}
					bChain.setAtom(connectedToRing, 0);
					bChain.setAtom(other1, 1);
					bChain.setOrder(b.getOrder());
//					if ((boolean) connectedToRing.getProperty(IS_LINKER_ATOM) == true)
//						bChain.setProperty(IS_LINKER, true);
//					else
//						bChain.setProperty(IS_SIDECHAIN, true);

					ac.addAtom(other1);

					bondToAdd.add(bChain);

					bondToRemove.add(b);
					counter++;
				}
			}
			else if (!b.getFlag(CDKConstants.ISINRING) && b.getAtom(1).getFlag(CDKConstants.ISINRING) 
					&& b.getAtom(0).getFlag(CDKConstants.ISINRING)) {
				if ((boolean) b.getAtom(0).getProperty(IS_LINKER_ATOM) && (boolean) b.getAtom(1).getProperty(IS_LINKER_ATOM)) {
					int counter2 = counter + 1;
					HAS_LINKER = true;
					
					IAtom inRing1 = b.getAtom(0);
					IAtom inRing2 = b.getAtom(1);
					
					inRing1.setProperty("connectorID", "R_"+counter);
					inRing2.setProperty("connectorID", "R_"+counter2);

					other1 = new Atom();
					other1.setSymbol("*");
					other1.setProperty("connection", inRing2);
					other1.setID("C_"+counter2);
					other1.setProperty("connectorID", "C_"+counter2);
					other1.setProperty(IS_LINKER_ATOM, true);
					other1.setProperty(IS_SIDECHAIN_ATOM, false);
					
					//bRing.setProperty("linkerConnector", true);


					other2 = new Atom();
					other2.setSymbol("*");
					other2.setProperty("connection", inRing1);
					other2.setID("C_"+counter);
					other2.setProperty("connectorID", "C_"+counter);
					other2.setProperty(IS_LINKER_ATOM, true);
					other1.setProperty(IS_SIDECHAIN_ATOM, false);

					IBond linker = new Bond();
					linker.setAtom(other1, 0);
					linker.setAtom(other2, 1);
					linker.setOrder(b.getOrder());

					ac.addAtom(other1);
					ac.addAtom(other2);

					bondToAdd.add(linker);

					bondToRemove.add(b);
					counter+=2;

				}
			}
		}
		
		//remove Bonds
		for (IBond b : bondToRemove)
			ac.removeBond(b);
		//add bonds
		for (IBond b : bondToAdd)
			ac.addBond(b);

		connections = new HashMap<IAtomContainer,List<IAtomContainer>>();

		//define an IAtomContainer for all independent fragments and return all these fragments in a atomContainer
		IAtomContainerSet set = FragmentUtils.makeAtomContainerSet(ac);

		//get index each atom (key:IAtom, value:index of the IATomContainer containing this atom)
		Map<IAtom,Integer> index = new HashMap<IAtom,Integer>();

		for (int i = 0; i < set.getAtomContainerCount(); i++) {
			for (IAtom a : set.getAtomContainer(i).atoms())
				index.put(a, i);
		}

		for (IAtomContainer ac2 : set.atomContainers()) {
//			try {
//				System.out.println(tools.makeSmiles(ac2, false, false));
//			} catch (CDKException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
			for (IAtom a : ac2.atoms()) {
				//it is a ring System
				if (a.getFlag(CDKConstants.ISINRING)) {
					ringSystems.add(ac2);
					break;
				}
				//find the corresponding ringSystem and add as connection
				if (a.getProperty("connection") != null) {
					IAtomContainer query = set.getAtomContainer(index.get(a.getProperty("connection")));

					if (connections.containsKey(query)) {
						List<IAtomContainer> list = connections.get(query);
						list.add(ac2);
						connections.put(query, list);
					}
					else {
						List<IAtomContainer> list = new ArrayList<IAtomContainer>();
						list.add(ac2);
						connections.put(query, list);
					}
//					try {
//						System.out.println(tools.makeSmiles(ac2, false, false));
//						System.out.println(a.getProperties());
//					} catch (CDKException e) {
//						// TODO Auto-generated catch block
//						e.printStackTrace();
//					}
					//a side chain has only one atom to connect
					if (a.getProperty(IS_SIDECHAIN_ATOM) != null) {
						if ((boolean)a.getProperty(IS_SIDECHAIN_ATOM) == true) {
							ac2.setProperty(IS_SIDECHAIN, true);
							break;
						}
					}	
					if (a.getProperty(IS_LINKER_ATOM) != null) {
						if ((boolean)a.getProperty(IS_LINKER_ATOM) == true)
							ac2.setProperty(IS_LINKER, true);
					}
				}
			}	
		}
	}

	
	/**
	 * generate a scaffold containing the ring system and and only the first connected bonds to the system <\br>
	 * the other atom in the group is replace by R group
	 * @param ac
	 * @return
	 */
	public IAtomContainer generateScaffold(IAtomContainer ac) {
		IAtomContainer scaffold = ac.getBuilder().newAtomContainer();

		List<IStereoElement> oldSe = getStereoElements(ac);
		Map<IChemObject,IStereoElement> focus = getFocus(oldSe);
		List<IStereoElement> seToAdd = new ArrayList<IStereoElement>();

		int counter = 1;

		for (IAtom a : ac.atoms()) {
			if (a.getFlag(CDKConstants.ISINRING) || (boolean) a.getProperty(IS_LINKER_ATOM)) {
				scaffold.addAtom(a);
				if (focus.containsKey(a))
					seToAdd.add(focus.get(a));
			}
		}
		for (IBond b : ac.bonds()) {
			if (b.getBegin().getFlag(CDKConstants.ISINRING) && b.getEnd().getFlag(CDKConstants.ISINRING)) {
				scaffold.addBond(b);
				if (focus.containsKey(b))
					seToAdd.add(focus.get(b));
			}
			else if (!b.getAtom(0).getFlag(CDKConstants.ISINRING) && b.getAtom(1).getFlag(CDKConstants.ISINRING)) {
				Atom a2 = new Atom();
				a2.setSymbol("R"+counter);
				counter++;
				a2.setID(b.getAtom(0).getID());
				b.setAtom(a2, 0);
				scaffold.addBond(b);
				if (focus.containsKey(b))
					seToAdd.add(focus.get(b));
			}
			else if (!b.getAtom(1).getFlag(CDKConstants.ISINRING)  && b.getAtom(0).getFlag(CDKConstants.ISINRING)) {
				Atom a2 = new Atom();
				a2.setSymbol("R"+counter);
				counter++;
				a2.setID(b.getAtom(1).getID());
				b.setAtom(a2, 1);
				scaffold.addBond(b);
				if (focus.containsKey(b))
					seToAdd.add(focus.get(b));
			}
		}

		scaffold.setStereoElements(seToAdd);

		return scaffold;
	}

	private Map<IChemObject,IStereoElement> getFocus(List<IStereoElement> seL) {
		Map<IChemObject,IStereoElement> focus = new HashMap<IChemObject,IStereoElement>();
		for (IStereoElement se : seL) {
			focus.put(se.getFocus(), se);
		}
		return focus;
	}

	private List<IStereoElement> getStereoElements(IAtomContainer ac) {
		List<IStereoElement> seL = new ArrayList<IStereoElement>();
		for (IStereoElement se : ac.stereoElements()) {
			seL.add(se);
		}
		return seL;
	}

	private void getOtherAtomsConnectedToARingSystem(Set<IAtom> sphere, IAtomContainer ac, Map<IAtom,IAtom> result) {

		Set<IAtom> newSphere = new HashSet<IAtom>();

		for (IAtom atom : sphere) {
			List<IBond> bonds = ac.getConnectedBondsList(atom);
			for (IBond bond : bonds) {
				IAtom nextAtom = bond.getConnectedAtom(atom);
				if (!nextAtom.getFlag(CDKConstants.ISINRING)) {
					if (!bond.getFlag(CDKConstants.VISITED)) {
						bond.setFlag(CDKConstants.VISITED, true);
					}
					if (!nextAtom.getFlag(CDKConstants.VISITED)) {
						if (!sphere.contains(nextAtom)) newSphere.add(nextAtom);
						nextAtom.setFlag(CDKConstants.VISITED, true);
					}
				}
				else
					result.put(atom, nextAtom);
			}
		}
		if (newSphere.size() > 0) {
			getOtherAtomsConnectedToARingSystem(newSphere, ac, result);
		}
	}

	public static List<IAtomContainer> createBricsFragment(IAtomContainer molecule, Map<IBond,String[]> bricsBonds,
			List<IBond> bricsBondsBroken) {
		List<IAtomContainer> result = new ArrayList<IAtomContainer>();
		while (bricsBondsBroken.size() < bricsBonds.size()) {
			if (bricsBondsBroken.size() == 0) {
				IBond cuttingBond = bricsBonds.keySet().iterator().next();
				addFragmentConnection(cuttingBond, molecule, result, bricsBonds.get(cuttingBond));
				bricsBondsBroken.add(cuttingBond);
			}
			else {
				for (IBond cuttingBond : bricsBonds.keySet()) {
					if (!bricsBondsBroken.contains(cuttingBond)) {
						for (IAtomContainer fragment : result) {
							if (fragment.contains(cuttingBond)) {
								addFragmentConnection(cuttingBond, fragment, result, bricsBonds.get(cuttingBond));
								bricsBondsBroken.add(cuttingBond);
								break;
							}
						}
					}
				}
			}
			//			System.out.println("t1 "+bricsFragments.size());
		}
		return result;
	}

	private static List<IAtomContainer> addFragmentConnection(IBond bond, IAtomContainer molecule, List<IAtomContainer> fragments, String[] def) {
		List<IAtomContainer> result = FragmentUtils.splitMolecule(molecule, bond);
		//put inside a function
		Atom a1 = new Atom("*");
		a1.setProperty("connectWithFragmentContainingAtom", bond.getAtom(1));
		a1.setProperty("def", def[1]);
		Bond b1 = new Bond();
		b1.setOrder(bond.getOrder());
		b1.setAtom(bond.getAtom(0), 0);
		b1.setAtom(a1, 1);
		result.get(0).addAtom(a1);
		result.get(0).addBond(b1);

		Atom a2 = new Atom("*");
		a2.setProperty("connectWithFragmentContainingAtom", bond.getAtom(0));
		a2.setProperty("def", def[0]);
		Bond b2 = new Bond();
		b2.setOrder(bond.getOrder());
		b2.setAtom(bond.getAtom(1), 0);
		b2.setAtom(a2, 1);
		result.get(1).addAtom(a2);
		result.get(1).addBond(b2);

		if (fragments.isEmpty())
			fragments.addAll(result);
		else {
			fragments.remove(molecule);
			fragments.addAll(result);
		}
		return fragments;
	}
	
	/**
	 * Make frameWork
	 * @param ac
	 * @return
	 */
	public static IAtomContainer getFramework(IAtomContainer ac) {
		IAtomContainer framework = DefaultChemObjectBuilder.getInstance().newAtomContainer();
		//make framework
		for (IBond b :ac.bonds()) {
			if (!b.getBegin().getSymbol().equals("R")) 
				framework.addAtom(b.getBegin());
			if (!b.getEnd().getSymbol().equals("R")) 
				framework.addAtom(b.getEnd());
			if (b.getProperty("connector") == null)
				framework.addBond(b);
		}
		return framework;
	}
	
	public boolean isHAS_LINKER() {
		return HAS_LINKER;
	}

	public void setHAS_LINKER(boolean hAS_LINKER) {
		HAS_LINKER = hAS_LINKER;
	}

	public List<IAtomContainer> getLinkers() {
		return linkers;
	}

	public List<IAtomContainer>  getSideChains() {
		return sideChains;
	}

	public List<IAtomContainer>  getRingSystems() {
		return ringSystems;
	}

	public Map<IAtomContainer, List<IAtomContainer>> getRingSystemLinkerConnections() {
		return ringSystemLinkerConnections;
	}

	public Map<IAtomContainer, List<IAtomContainer>> getRingSystemSideChainConnections() {
		return ringSystemSideChainConnections;
	}
	
	public Map<IAtomContainer, List<IAtomContainer>> getConnections() {
		return connections;
	}

	public Map<String, String> getConnectedAtoms() {
		return connectedAtoms;
	}

	public Map<String, IBond> getRingToConnector() {
		return ringToConnector;
	}
	
}
