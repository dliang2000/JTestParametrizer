package ca.uwaterloo.jrefactoring.match;

import java.util.*;

public class InputMethodsGroup {
	private int inputMethodsSize;
	private Set<InputMethodInstance> inputMethodInstances = new LinkedHashSet<InputMethodInstance>();
	
	public int getInputMethodsSize() {
		return inputMethodsSize;
	}

	public void setInputMethodsSize(int inputMethodsSize) {
		this.inputMethodsSize = inputMethodsSize;
	}
	
	public void addInputMethod(InputMethodInstance inputMethodInstance) {
		inputMethodInstances.add(inputMethodInstance);
	}
	
	public List<InputMethodInstance> getInputMethodInstances() {
		List<InputMethodInstance> copyToReturn = new ArrayList<InputMethodInstance>();
		for (InputMethodInstance inputMethodInstance : inputMethodInstances)
			copyToReturn.add(inputMethodInstance);
		return copyToReturn;
	}
}
