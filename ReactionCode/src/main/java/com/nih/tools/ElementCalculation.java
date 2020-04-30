package com.nih.tools;

public class ElementCalculation {

	public static int calculateValence(String symbol) {
		Element e = null;
		
		//if (symbol.equals("C")) e = Element.C;
		//int valence = 4;

		e = Element.valueOfIgnoreCase(symbol);
		//System.out.println(e.toString()+" "+e.getCommonValence() + " " +e.getMinimumValence()+" "+e.getMaximumValence());
		
        return e.getCommonValence();
	}
	
	public static int calculateMass(String symbol) {
		Element e = null;
		
		//if (symbol.equals("C")) e = Element.C;
		//int valence = 4;

		e = Element.valueOfIgnoreCase(symbol);
		//System.out.println(e.toString()+" "+e.getCommonValence() + " " +e.getMinimumValence()+" "+e.getMaximumValence());
		
        return (int) e.getAtomicMass();
	}
	
	public static boolean isHalogen(String symbol) {
		boolean isHalogen = false;
		Element e = null;
		
		//if (symbol.equals("C")) e = Element.C;
		//int valence = 4;

		e = Element.valueOfIgnoreCase(symbol);
		if (e.getElementType() == ElementType.HALOGEN) isHalogen = true;
		
		//System.out.println(e.toString()+" "+e.getCommonValence() + " " +e.getMinimumValence()+" "+e.getMaximumValence());
		
        return isHalogen;
	}
}
