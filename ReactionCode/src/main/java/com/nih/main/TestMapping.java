package com.nih.main;

import java.io.IOException;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IReaction;

import com.nih.mapping.Mapping;

public class TestMapping {
	
	

	public static void main(String[] args) throws IOException, CDKException, InterruptedException {
		// TODO Auto-generated method stub
//		String smirks = "FC(F)(F)C1=CC=C(Cl)N=C1.OC[@H](C)C2=CC(Br)=CC=C2.[H-]>>OC[@@H](CC3=CC=C(C=N3)C(F)(F)F)C4=CC(Br)=CC=C4.[Cl-].[HH]";
//		String smirks = "FC(F)(F)C1=CC=C(Cl)N=C1.O=C(C)C2=CC(Br)=CC=C2>>O=C(CC3=CC=C(C=N3)C(F)(F)F)C4=CC(Br)=CC=C4.HCl";
		String smirks = "FC(F)(F)C1=CC=C(Cl)N=C1.O=C(C)C2=CC(Br)=CC=C2.[H-]>>O=C(CC3=CC=C(C=N3)C(F)(F)F)C4=CC(Br)=CC=C4.[Cl-].[HH]";
//		String smirks = "FC(F)(F)C1=CC=C(Cl)N=C1.O=C(C)C2=CC(Br)=CC=C2.[H-]>>O=C(CC3=CC=C(C=N3)C(F)(F)F)C4=CC(Br)=CC=C4.[Cl-].HH";
		
		Mapping mapping = new Mapping();
		IReaction reaction = mapping.mapReaction(smirks);
		System.out.println(reaction);
		
	}
	
	
}
