package com.nih.fragments;

/*Adapted from BRICS in RDKIT*/

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.lang3.ArrayUtils;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smarts.SmartsPattern;

public class BRICS {

	Map<String, String> environs = new HashMap<String,String>();
	List<String[][]> reactionDefs = new ArrayList<String[][]>();
	List<String[]> bondMatchers = new ArrayList<String[]>();

	public void initialize() {
		environs.put("L1","[C;D3]([#0,#6,#7,#8])(=O)");
		  /* After some discussion, the L2 definitions ("N.pl3" in the original
		   paper) have been removed and incorporated into a (almost) general
		   purpose amine definition in L5 ("N.sp3" in the paper).
		  
		   The problem is one of consistency.
		      Based on the original definitions you should get the following
		      fragmentations:
		        C1CCCCC1NC(=O)C -> C1CCCCC1N[2*].[1*]C(=O)C
		        c1ccccc1NC(=O)C -> c1ccccc1[16*].[2*]N[2*].[1*]C(=O)C
		      This difference just didn't make sense to us. By switching to
		      the unified definition we end up with:
		        C1CCCCC1NC(=O)C -> C1CCCCC1[15*].[5*]N[5*].[1*]C(=O)C
		        c1ccccc1NC(=O)C -> c1ccccc1[16*].[5*]N[5*].[1*]C(=O)C
		        */
		environs.put("L2", "[N;!R;!D1;!$(N=*)]-;!@[#0,#6]");
//		this one turned out to be too tricky to define above, so we set it off in its own definition:
//		environs.put("L2a", "[N;D3;R;$(N(@[C;!$(C=*)])@[C;!$(C=*)])]");
		environs.put("L3", "[O;D2]-;!@[#0,#6,#1]");
		environs.put("L4", "[C;!D1;!$(C=*)]-;!@[#6]");
//		environs.put("L5", "[N;!D1;!$(N*!-*);!$(N=*);!$(N-[!C;!#0])]-[#0,C]");
		environs.put("L5", "[N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]");
		environs.put("L6", "[C;D3;!R](=O)-;!@[#0,#6,#7,#8]");
		environs.put("L7a", "[C;D2,D3]-[#6]");
		environs.put("L7b", "[C;D2,D3]-[#6]");
//		environs.put("L8", "[C;!R;!D1]-;!@[#6]");
		environs.put("L8", "[C;!R;!D1;!$(C!-*)]");
		environs.put("L9", "[n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]");
		environs.put("L10", "[N;R;$(N(@C(=O))@[C,N,O,S])]");
		environs.put("L11", "[S;D2](-;!@[#0,#6])");
		environs.put("L12", "[S;D4]([#6,#0])(=O)(=O)");
		environs.put("L13", "[C;$(C(-;@[C,N,O,S])-;@[N,O,S])]");
		environs.put("L14", "[c;$(c(:[c,n,o,s]):[n,o,s])]");
		environs.put("L14b", "[c;$(c(:[c,n,o,s]):[n,o,s])]");
		environs.put("L15", "[C;$(C(-;@C)-;@C)]");
		environs.put("L16", "[c;$(c(:c):c)]");
		environs.put("L16b", "[c;$(c(:c):c)]");

		 String L1[][] = { 
				 {"1","3","-"},
				 {"1","5","-"},
				 {"1","10","-"} 
				 };
		 String L3[][] = { 
				 {"3","4","-"},
				 {"3","13","-"},
				 {"3","14","-"},
				 {"3","15","-"},
				 {"3","16","-"} 
				 };
		 String L4[][] = { 
				 {"4","5","-"},
				 {"4","11","-"} 
				 };
		 String L5[][] = { 
				 {"5","12","-"},
				 {"5","14","-"},
				 {"5","16","-"},
				 {"5","13","-"},
				 {"5","15","-"} 
				 };
		 String L6[][] = { 
				 {"6","13","-"},
				 {"6","14","-"},
				 {"6","15","-"},
				 {"6","16","-"} 
				 };
		 String L7[][] = { 
				 {"7a","7b","="}
				 };
		 String L8[][] = { 
				 {"8","9","-"},
				 {"8","10","-"},
				 {"8","13","-"},
				 {"8","14","-"},
				 {"8","15","-"},
				 {"8","16","-"} 
				 };
		 String L9[][] = { 
				 {"9","13","-"}, //not in original paper
				 {"9","14","-"}, //not in original paper
				 {"9","15","-"},
				 {"9","16","-"} 
				 };
		 String L10[][] = { 
				 {"10","13","-"}, 
				 {"10","14","-"}, 
				 {"10","15","-"},
				 {"10","16","-"} 
				 };
		 String L11[][] = { 
				 {"11","13","-"}, 
				 {"11","14","-"}, 
				 {"11","15","-"},
				 {"11","16","-"} 
				 };
//		 String L12[][] = { 
//				 {"","",""}, //none
//				 };
		 String L13[][] = { 
				 {"13","14","-"}, 
				 {"13","15","-"}, 
				 {"13","16","-"} 
				 };
		 String L14[][] = { 
				 {"14","14","-"}, //not in original paper
				 {"14","15","-"}, 
				 {"14","16","-"} 
				 };
		 String L15[][] = { 
				 {"15","16","-"} 
				 };
		 String L16[][] = { 
				 {"16","16","-"} //not in original paper
				 };
		 reactionDefs.add(L1);
		 reactionDefs.add(L3);
		 reactionDefs.add(L4);
		 reactionDefs.add(L5);
		 reactionDefs.add(L6);
		 reactionDefs.add(L7);
		 reactionDefs.add(L8);
		 reactionDefs.add(L9);
		 reactionDefs.add(L10);
		 reactionDefs.add(L11);
		 reactionDefs.add(L13);
		 reactionDefs.add(L14);
		 reactionDefs.add(L15);
		 reactionDefs.add(L16);
		 
		 for (String[][] def : reactionDefs) {
			 for (int i = 0; i < def.length; i++) {
				 String e1 = environs.get("L"+def[i][0]);
				 String e2 = environs.get("L"+def[i][1]);
				 String bType = def[i][2];
				 String pattern = "[$(" + e1 +")]" + bType + ";!@[$(" + e2 + ")]";
				 String matcher[] = {def[i][0],def[i][1],bType,pattern};
				 bondMatchers.add(matcher);
			 } 
		 }
	}
	
	public Map<IBond,String[]> FindBRICSBonds(IAtomContainer molecule) {
		Set<IBond> bDone = new HashSet<IBond>();
		Map<IBond,String[]> result = new HashMap<IBond,String[]>();
		Map<String,Mappings> envMatches = new HashMap<String,Mappings>();
		for (Entry<String,String> e: environs.entrySet()) {
			 Pattern ptrn = SmartsPattern.create(e.getValue());
			 envMatches.put(e.getKey(), ptrn.matchAll(molecule));
		}
		for (String[] bMatcher : bondMatchers) {
			if (!envMatches.containsKey("L"+bMatcher[0]) || !envMatches.containsKey("L"+bMatcher[1]))
				continue;
			Pattern ptrn = SmartsPattern.create(bMatcher[3]);
			int atom2Index = bMatcher[3].indexOf(";!@[$(")+1;
//			System.out.println(bMatcher[3].substring(atom2Index));
			Pattern ptrn2 = SmartsPattern.create(bMatcher[3].substring(atom2Index));
			Iterable<IAtomContainer> iacs = ptrn2.matchAll(molecule).toSubstructures();
			for (IAtomContainer ac : ptrn.matchAll(molecule).toSubstructures()) {
				for (IBond bMatch : ac.bonds()) {
					if (!bDone.contains(bMatch)) {
//						System.out.println(bMatch);
						//add fragment connection information throw atom of the target fragment
						bMatch.getBegin().setProperty("connectWithFragmentContainingAtom", bMatch.getEnd());
						bMatch.getEnd().setProperty("connectWithFragmentContainingAtom", bMatch.getBegin());
						//check if second smart pattern match with 2d atom in the matched bond
						for (IAtomContainer iac : iacs) {
							for (IAtom a : iac.atoms()) {
								//reverse def array if 2d smart pattern match with 1st atom in the bond
								if (bMatch.getAtom(0).equals(a)) {
									String[] arr = {bMatcher[1], bMatcher[0]};
									result.put(bMatch,arr);
								}
								//add def array corresponding to def[0] -> atom 0 in Bond and def[1] -> atom 1 in Bond
								if (bMatch.getAtom(1).equals(a))
									result.put(bMatch,Arrays.copyOfRange(bMatcher, 0, 2));
							}
						}
//						result.put(bMatch,Arrays.copyOfRange(bMatcher, 0, 2));
					}
				}
			}
		}
		return result;
	}
}
