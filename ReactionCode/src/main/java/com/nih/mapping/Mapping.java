package com.nih.mapping;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomContainerSetManipulator;

public class Mapping {
	int counter = 0;
	int counter2 = 0;
	int atomsInReactant = 0;
	String newsmirks = "";
	List<IAtom> index = new ArrayList<IAtom>();
	IReaction reaction;
	List<String> reactionCenter = new ArrayList<String>();
	Map<IAtom,IAtomContainer> index2;
	Map<IAtom,IAtomContainer> explicitHydrogensInReactants;
	Map<IAtom,IAtomContainer> explicitHydrogensInProducts;
	

	public IReaction mapReaction(String smirks) throws IOException, CDKException, InterruptedException {
		IReaction reaction = parseReactionSmiles(smirks);
		return mapReaction(reaction);
	}
	
	public IReaction mapReaction(IReaction reac) throws IOException, CDKException, InterruptedException {
		reaction = reac;
		//prepare reactants and products by storing the non H atom in index as a function of their index and make new SMIRKS
		explicitHydrogensInReactants = prepareAtomContainer(reaction.getReactants());
		newsmirks += ">>";
		explicitHydrogensInProducts = prepareAtomContainer(reaction.getProducts());

//		System.out.println(newsmirks);
		
		//config 
		String cmd = "/home/delanneev/Downloads/AtomMap-master/map_reaction";
		String directory = "/home/delanneev/Downloads/AtomMap-master/";

		//make Mapping
		InputStream result = jaworskiMapping(cmd, directory);

		//apply it to CDK Object
		boolean status = applyMappingToCDKReaction(result);
		
		if (status == false)
			return null;
		else
			return reaction;
	}
	
	private IReaction parseReactionSmiles(String smirks) {
		IReaction reaction = null;
		try {
			SmilesParser   sp  = new SmilesParser(SilentChemObjectBuilder.getInstance());
			reaction   = sp.parseReactionSmiles(smirks);
		} catch (InvalidSmilesException e) {
			System.err.println(e.getMessage());
		}
		return reaction;
	}
	
	private Map<IAtom,IAtomContainer> prepareAtomContainer(IAtomContainerSet set) throws CDKException {
		index2 = new HashMap<IAtom,IAtomContainer>();
		Map<IAtom,IAtomContainer> explicitHydrogens = new HashMap<IAtom,IAtomContainer>();
		SmilesGenerator smigen = new SmilesGenerator(SmiFlavor.Isomeric);
		for (IAtomContainer ac : set.atomContainers()) {
//			AtomContainerSetManipulator.getTotalHydrogenCount(set);
//			AtomContainerManipulator.convertImplicitToExplicitHydrogens(ac);
			if (newsmirks.length() > 0) {
				newsmirks += ".";
			}
			int[] order = new int[ac.getAtomCount()];
			String smi = smigen.create(ac, order);
			newsmirks += smi;

//			System.out.println(smi + ac.getAtomCount());
			for (int o : order) {
				IAtom atom = ac.getAtom(o);
				counter2++;
				atom.setID(counter2+"");
				if (!atom.getSymbol().equals("H")){
					index.add(atom);
					index2.put(atom, ac);
					counter++;
				}
				else {
					index2.put(atom, ac);
					explicitHydrogens.put(atom, ac);
				}
			}
		}
		return explicitHydrogens;
	}

	private InputStream jaworskiMapping(String cmd, String directory) throws InterruptedException, IOException {
		ProcessBuilder builder = new ProcessBuilder();
		builder.directory(new File(directory));
		builder.command(cmd, newsmirks);
		Process process = builder.start();
		process.waitFor();
		
		int exitValue = process.exitValue();
//		System.out.println("exitValue = " + exitValue);
		
		InputStream stdout = process.getInputStream();
//		InputStream stderr = process.getErrorStream();
//		
//		BufferedReader errorred = new BufferedReader(new InputStreamReader(stderr));
//		String line2 = null;
//		while ((line2 = errorred.readLine()) != null) {
//			System.out.println(line2);
//		}
		
		return stdout;
	}
	
	private boolean applyMappingToCDKReaction(InputStream stream) throws IOException {
		BufferedReader reader = new BufferedReader(new InputStreamReader(stream));
		String line1 = null;
		int atomCounter = -1;
		atomsInReactant = 0;
		
		while ((line1 = reader.readLine()) != null) {
			Map<String,String> indexAtom = new HashMap<String,String>();
			Map<Integer,IAtom> mappingOfSmilesWithID = new TreeMap<Integer,IAtom>();
			String reactionSmiles = null;
			if (line1.contains("<reaction>")) {
				do {
					line1 = reader.readLine();
					if (line1.contains("<smiles>")) {
						reactionSmiles = line1.replace("<smiles>", "")
								.replace("</smiles>", "")
								.replaceAll("\\s+", "")
								.replaceAll("&gt;", ">");
					}
					if (line1.contains("<messages>")) {
						System.err.println("SMIRKS" + reactionSmiles);
						System.err.println(" \t" +reader.readLine().replace("<message>", "")
								.replace("</message>", "").
								replaceAll("\\s+", "").replaceAll("&quot;", ""));
						return false;
					}
					
				} while (!line1.contains("<solution>"));
				
				do {
					line1 = reader.readLine();
					String reactantsSMILESWithID = line1.replace("<reactants>", "").replace("</reactants>", "").replaceAll("\\s+", "");
					String productsSMILESWithID = line1.replace("<products>", "").replace("</products>", "").replaceAll("\\s+", "");
	
					indexAtom.putAll(makeAtomIndex(reactantsSMILESWithID));
					atomsInReactant = indexAtom.size();
					indexAtom.putAll(makeAtomIndex(productsSMILESWithID));
					
//					System.out.println(reactantsSMILESWithID);
//					System.out.println(productsSMILESWithID);
	//				System.out.println(indexAtom);
					
					if (line1.contains("<atoms>")) {
						String atID = null;
						String visible;
						do {
							line1 = reader.readLine();
							
							if (line1.contains("atom id=")) {
								atID = line1.split("\"")[1];
								visible = line1.split("\"")[3];
							}
							if (line1.contains("<label>")) {
								if (!indexAtom.get(atID).equals("H")) {
									atomCounter += 1;
									String label = line1.replace("<label>", "").replace("</label>", "").replaceAll("\\s+", "");
									IAtom atom = index.get(atomCounter);
									atom.setProperty(CDKConstants.ATOM_ATOM_MAPPING, Integer.parseInt(label));
									mappingOfSmilesWithID.put(Integer.parseInt(atID), atom);
								}
								
							}
						} while (!line1.contains("</atoms>"));
					}
					if (line1.contains("<bonds>")) {
						String at1ID = null;
						String at2ID = null;

						do {
							line1 = reader.readLine();
							if (line1.contains("id1=")) {
								at1ID = line1.split("\"")[1].replaceAll("\\s+", "");
								at2ID = line1.split("\"")[3].replaceAll("\\s+", "");
								reactionCenter.add(at1ID);
								reactionCenter.add(at2ID);
								reactionCenter.add(at1ID+"-"+at2ID);
								reactionCenter.add(at2ID+"-"+at1ID);
							}
							
						} while (!line1.contains("</bonds>"));
					}
				} while (!line1.contains("</solution>"));
				mapExplicitHydrogen(indexAtom, mappingOfSmilesWithID);	
			}
		}
		reader.close();
		
		
		
//		for (IAtom a : index) {
//			System.out.println(a.getSymbol() + " " +a.getProperties());
//		}
		
		
		return true;
	}
	
	private void mapExplicitHydrogen(Map<String,String> indexAtom, Map<Integer,IAtom> mappingOfSmilesWithID) {
		List<IAtom> processed = new ArrayList<IAtom>();
		for (Entry<IAtom,IAtomContainer> e2 : explicitHydrogensInProducts.entrySet()) {
			IAtom hydrogen2 = e2.getKey();
			if (processed.contains(hydrogen2)) {
				continue;
			}
			IAtomContainer molecule2 = e2.getValue();
			List<IAtom> con = molecule2.getConnectedAtomsList(hydrogen2);
			if (con.size() == 0) {
				if (hydrogen2.getImplicitHydrogenCount() != null) {
					if (hydrogen2.getImplicitHydrogenCount() == 1) {
						for (Entry<IAtom,IAtomContainer> e1 : explicitHydrogensInReactants.entrySet()) {
							IAtom hydrogen1 = e1.getKey();
							if (processed.contains(hydrogen1)) {
								continue;
							}
							IAtomContainer molecule1 = e1.getValue();
							List<Integer> con1 = new ArrayList<Integer>();
							molecule1.getConnectedAtomsList(hydrogen1)
									.forEach(a -> con1.add(a.getProperty(CDKConstants.ATOM_ATOM_MAPPING)));

							if (con1.size() == 0) {
								int aam = Integer.parseInt(hydrogen1.getID());
								hydrogen2.setProperty(CDKConstants.ATOM_ATOM_MAPPING, aam);
								hydrogen1.setProperty(CDKConstants.ATOM_ATOM_MAPPING, aam);
								
								processed.add(hydrogen2);
								processed.add(hydrogen1);

								break;
							}
						}
					}
				}
			}
			else {
				IAtom other = con.get(0);
				for (Entry<IAtom,IAtomContainer> e1 : explicitHydrogensInReactants.entrySet()) {
					IAtom hydrogen1 = e1.getKey();
					if (processed.contains(hydrogen1)) {
						continue;
					}
					IAtomContainer molecule1 = e1.getValue();
					List<Integer> con1 = new ArrayList<Integer>();
					molecule1.getConnectedAtomsList(hydrogen1)
							.forEach(a -> con1.add(a.getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
					
					if (con1.size() > 0) {
						if (other.getProperty(CDKConstants.ATOM_ATOM_MAPPING) == con1.get(0)) {
							int aam = Integer.parseInt(hydrogen1.getID());
							hydrogen2.setProperty(CDKConstants.ATOM_ATOM_MAPPING, aam);
							hydrogen1.setProperty(CDKConstants.ATOM_ATOM_MAPPING, aam);
							processed.add(hydrogen2);
							processed.add(hydrogen1);
							break;
						}
					}
				}
				if (!processed.contains(hydrogen2) && other.getSymbol().equals("H")) {
					for (Entry<IAtom,IAtomContainer> e1 : explicitHydrogensInReactants.entrySet()) {
						IAtom hydrogen1 = e1.getKey();
						if (processed.contains(hydrogen1)) {
							continue;
						}
						IAtomContainer molecule1 = e1.getValue();
						List<Integer> con1 = new ArrayList<Integer>();
						molecule1.getConnectedAtomsList(hydrogen1)
								.forEach(a -> con1.add(a.getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
						if (con1.size() == 0) {
							int aam = Integer.parseInt(hydrogen1.getID());
							hydrogen2.setProperty(CDKConstants.ATOM_ATOM_MAPPING, aam);
							hydrogen1.setProperty(CDKConstants.ATOM_ATOM_MAPPING, aam);
							processed.add(other);
							processed.add(hydrogen2);
							processed.add(hydrogen1);
							
							molecule2.removeBond(molecule2.getConnectedBondsList(hydrogen2).get(0));
							molecule2.removeAtom(other);
							hydrogen2.setImplicitHydrogenCount(1);
							System.out.println("Test");
							break;
						}
					}
				}
			}
		}
		

		
//		for (String rcID1 : reactionCenter) {
//			if (rcID1.contains("-")) {
//				continue;
//			}
//			String symbol = indexAtom.get(rcID1);
//			if (symbol.equals("H") && Integer.parseInt(rcID1) > atomsInReactant) {
//				String rcID2 = null;
//				for (String str : reactionCenter) {
//					if (str.contains("rcID1-")) {
//						rcID2 = str.replace("rcID1-", "");
//					}
//				}
//				//look for HH
//				if (indexAtom.get(rcID2).equals("H")) {
//					for (Entry<IAtom,IAtomContainer> e2 : explicitHydrogensInProducts.entrySet()) {
//						IAtom hydrogen2 = e2.getKey();
//						if (processed.contains(hydrogen2)) {
//							continue;
//						}
//						IAtomContainer molecule2 = e2.getValue();
//						List<IAtom> con = molecule2.getConnectedAtomsList(hydrogen2);
//						if (con.size() == 0 && hydrogen2.getImplicitHydrogenCount() == 1) {
//							
//						}
//
//					}
//				}
//				int connectedHeteroAtomInt = Integer.parseInt(rcID1);
//				while (mappingOfSmilesWithID.containsKey(connectedHeteroAtomInt)) {
//					connectedHeteroAtomInt--;
//				}
//				IAtom connectedHeteroAtom = mappingOfSmilesWithID.get(connectedHeteroAtomInt);
//				List<IAtom> hydrogens = new ArrayList<IAtom>();
//				List<IAtom> con = index2.get(connectedHeteroAtom).getConnectedAtomsList(connectedHeteroAtom);
//			}
//			//an hydrogen bond is broken
//			if (symbol.equals("H") && Integer.parseInt(rcID1) <= atomsInReactant) {
//				int connectedHeteroAtomInt = Integer.parseInt(rcID1);
//				while (mappingOfSmilesWithID.containsKey(connectedHeteroAtomInt)) {
//					connectedHeteroAtomInt--;
//				}
//				IAtom connectedHeteroAtom = mappingOfSmilesWithID.get(connectedHeteroAtomInt);
//				List<IAtom> hydrogens = new ArrayList<IAtom>();
//				List<IAtom> con = index2.get(connectedHeteroAtom).getConnectedAtomsList(connectedHeteroAtom);
//				for (IAtom a : con) {
//					if ( a.getSymbol().equals("H")) {
//						hydrogens.add(a);
//					}
//				}
//				if (hydrogens.isEmpty()) {
//					
//				}
//				//look for a non connected atom
//				for (Entry<IAtom,IAtomContainer> e2 : explicitHydrogensInProducts.entrySet()) {
//					IAtom hydrogen2 = e2.getKey();
//					if (processed.contains(hydrogen2)) {
//						continue;
//					}
//					IAtomContainer molecule2 = e2.getValue();
//				}
//			}
//		}
//		
//		counter++;
//		List<IAtom> processed = new ArrayList<IAtom>();
//		System.out.println(explicitHydrogensInReactants.size());
//		System.out.println(explicitHydrogensInProducts.size());
//		for (Entry<IAtom,IAtomContainer> e1 : explicitHydrogensInReactants.entrySet()) {
//			IAtom hydrogen1 = e1.getKey();
//			IAtomContainer molecule1 = e1.getValue();
//			List<Integer> con1 = new ArrayList<Integer>();
//			molecule1.getConnectedAtomsList(hydrogen1)
//					.forEach(a -> con1.add(a.getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
//			
//			//look for an alone hydrogen or an hydrogen involved in a new formed bond
//			if (con1.size() == 0) {
//				
//			}
//			else {
//				for (Entry<IAtom,IAtomContainer> e2 : explicitHydrogensInProducts.entrySet()) {
//					IAtom hydrogen2 = e2.getKey();
//					if (processed.contains(hydrogen2)) {
//						continue;
//					}
//					IAtomContainer molecule2 = e2.getValue();
//					List<Integer> con2 = new ArrayList<Integer>();
//					molecule2.getConnectedAtomsList(hydrogen2)
//							.forEach(a -> con2.add(a.getProperty(CDKConstants.ATOM_ATOM_MAPPING)));
//					if (molecule1.getConnectedAtomsList(hydrogen1).size() > 0 &&
//					molecule2.getConnectedAtomsList(hydrogen2).size() > 0) {
//						System.out.println(molecule1.getConnectedAtomsList(hydrogen1).get(0).getProperties());
//						System.out.println(molecule2.getConnectedAtomsList(hydrogen2).get(0).getProperties());
//					}
//					
//					System.out.println(hydrogen1.getID() + " " + con1 + " " +hydrogen2.getID() + " " +  con2);
//					if (!con1.isEmpty() && con1.containsAll(con2)) {
////						System.out.println(con1 + " " + con2);
//						hydrogen1.setProperty(CDKConstants.ATOM_ATOM_MAPPING, counter);
//						hydrogen2.setProperty(CDKConstants.ATOM_ATOM_MAPPING, counter);
//						processed.add(hydrogen1);
//						processed.add(hydrogen2);
//						counter++;
//						break;
//					}
//				}
//			}
//			
//		}
		
		
		//To check if all explicit hydrogen have been mapped
//		for (Entry<IAtom,IAtomContainer> e2 : explicitHydrogensInReactants.entrySet()) {
//			IAtom hydrogen2 = e2.getKey();
//			if (!processed.contains(hydrogen2)) {
//				System.out.println("TODO " +e2.getKey());
//			}
//		}
//		for (Entry<IAtom,IAtomContainer> e2 : explicitHydrogensInProducts.entrySet()) {
//			IAtom hydrogen2 = e2.getKey();
//			if (!processed.contains(hydrogen2)) {
//				System.out.println("TODO2 " +e2.getKey());
//			}
//		}
	}
	
	private static Map<String,String> makeAtomIndex(String str) {
		Map<String,String> result = new HashMap<String,String>();
		char[] chars = str.toCharArray();
		for (int i = 0; i < chars.length; i++) {
			Character c = chars[i];
			if (c == '[') {
				String  symbol = "";
				String  id = "";
				i++;
				do {
					symbol += chars[i];
					i++;
					c = chars[i];
				} while (c != ':');
				i++;
				do {
					id += chars[i];
					i++;
					c = chars[i];
				} while (c != ']');
				result.put(id, symbol);
			}
		}
		return result;
	}
}
