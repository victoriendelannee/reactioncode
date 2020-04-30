package com.nih.main;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IReaction;

import com.nih.codes.DecodeReactionCode;
import com.nih.tools.tools;
import com.nih.transformer.Transformer;

public class TestTransformer {
	public static void main(String[] args) throws CDKException, IOException, CloneNotSupportedException {
		DecodeReactionCode decoder = new DecodeReactionCode();
		Transformer transformer = new Transformer();
		IReaction ptrn = decoder.decode("0:908()[1]507()[1]506(21GH01GG)[1]|");
		String smi = "[C:21]([CH3:22])([CH3:23])([CH3:24])[N:25]=[C:26]=[O:27].[H-:19].[Na+:20].[O:28]1[CH2:29][CH2:30][CH2:31][CH2:32]1.[OH:1][c:2]1[cH:3][cH:4][c:5]([CH2:8][CH2:9][O:10][c:11]2[cH:12][cH:13][c:14]([CH2:15][OH:16])[cH:17][cH:18]2)[cH:6][cH:7]1";
		transformer.setCheckValence(true);
		transformer.setStoichiometry(1);
		List<IReaction> trans = transformer.applyTranform2(smi, ptrn);
		System.out.println(tools.makeSmiles(ptrn, false));
		//System.out.println(tools.makeSmiles(transformer.applyTranform("[NH2:1][CH:6]1[CH2:11][CH2:17]([NH2:26])[CH2:22][CH2:18][CH2:12]1.[Cl:2][cH:3]1[n:4]([nH:7][cH:8]2[cH:5]1[cH2:10][cH2:16][cH:19]([cH2:13]2)[F:23])-[cH:9]3[cH2:14][cH2:20][cH:24]([cH2:21][cH2:15]3)[Cl:25]", ptrn), false));
		System.out.println(trans.size());
		for (IReaction reaction : trans) {
			System.out.println(tools.makeSmiles(reaction, true));
		}
		List<String> reactionCodes = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new FileReader("/home/delanneev/Downloads/zzz_reactionCodes_converted.txt")); 
		String line;
		while ((line = br.readLine()) != null) {
			reactionCodes.add(line);
		}
		br = new BufferedReader(new FileReader("/home/delanneev/Downloads/zzz_reagents.txt")); 
		int cpt = 0;
		while ((line = br.readLine()) != null) {
			System.out.println("------------------");
			System.out.println(reactionCodes.get(cpt));
			IReaction ptrn2 = decoder.decode(reactionCodes.get(cpt));
			System.out.println(line);
			System.out.println(tools.makeSmiles(ptrn2, false));
			List<IReaction> trans2 = transformer.applyTranform2(line, ptrn2);
			System.out.println(trans2.size());
			cpt++;

		}
		}
}
