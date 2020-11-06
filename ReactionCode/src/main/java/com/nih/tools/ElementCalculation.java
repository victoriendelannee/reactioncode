package com.nih.tools;

public class ElementCalculation {

	public static int calculateValence(Element elem) {
        return elem.getCommonValence();
	}
	
	public static int calculateValence(String symbol) {
		Element e = null;
		
		e = Element.valueOfIgnoreCase(symbol);
		
        return e.getCommonValence();
	}
	
	
	public static int calculateMinimumValence(String symbol) {
		Element e = null;

		e = Element.valueOfIgnoreCase(symbol);		
        return e.getMinimumValence();
	}
	
	public static int calculateMaximumValence(String symbol) {
		Element e = null;

		e = Element.valueOfIgnoreCase(symbol);
		
        return e.getMaximumValence();
	}
	
	public static int calculateMass(String symbol) {
		Element e = null;

		e = Element.valueOfIgnoreCase(symbol);
		
        return (int) e.getAtomicMass();
	}
	
	public static boolean isHalogen(String symbol) {
		boolean isHalogen = false;
		Element e = null;

		e = Element.valueOfIgnoreCase(symbol);
		if (e.getElementType() == ElementType.HALOGEN) isHalogen = true;
				
        return isHalogen;
	}
}
