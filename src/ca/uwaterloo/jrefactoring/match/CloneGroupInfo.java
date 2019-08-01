package ca.uwaterloo.jrefactoring.match;

import org.eclipse.core.filebuffers.FileBuffers;
import org.eclipse.core.filebuffers.ITextFileBuffer;
import org.eclipse.core.filebuffers.ITextFileBufferManager;
import org.eclipse.core.filebuffers.LocationKind;
import org.eclipse.core.runtime.CoreException;
import org.eclipse.core.runtime.IPath;
import org.eclipse.jdt.core.*;
import org.eclipse.jface.text.BadLocationException;
import org.eclipse.jface.text.IDocument;

import java.util.*;

public class CloneGroupInfo {
	
	public enum AnalysisStatus {
		/** Happens when 
		 * <ul>
		 * <li>At least one of the ASTs didn't have any nodes</li>
		 * <li>SystemObject.getMethodObject() cannot find either iMethod1 or iMethod2</li>,
		 * <li>The getMethodBody() of either methodObject1 or methodObject2 (resulting from calling SystemObject.getMethodObject()) is null</li>
		 * </ul>
		 */
		NOT_ANALYZED,
		/** The bottom-up subtree mapping didn't find any common nesting structure */
		NO_COMMON_SUBTREE_FOUND,
		/** Normal */
		NORMAL


	}
	
	private String projectName;
	private int cloneGroupID;
	private int cloneGroupSize;
	private Set<CloneItem> CloneItems = new LinkedHashSet<CloneItem>();
	
	private List<Long> subtreeMatchingTimeList = new ArrayList<>();
	private List<Integer> numberOfStatementsToBeRefactored = new ArrayList<>();
	private List<Integer> numberOfNodeComparisons = new ArrayList<>();
	private AnalysisStatus status;
	private List<List<PDGSubTreeMapperInfo>> subtreeMapperListList = new ArrayList<List<PDGSubTreeMapperInfo>>();
	private Set<String> testPackages;
	private Set<String> testSourceFolders;

	public String getProjectName() {
		return projectName;
	}

	public void setProjectName(String projectName) {
		this.projectName = projectName;
	}
	
	public int getCloneGroupID() {
		return cloneGroupID;
	}

	public void setCloneGroupID(int cloneGroupID) {
		this.cloneGroupID = cloneGroupID;
	}
	
	public int getCloneGroupSize() {
		return cloneGroupSize;
	}

	public void setCloneGroupSize(int cloneGroupSize) {
		this.cloneGroupSize = cloneGroupSize;
	}
	
	public void addClone(CloneItem cloneItem) {
		if (!cloneItem.isClassLevelClone()) {
			CloneItems.add(cloneItem);
			cloneItem.setCloneGroupBelonging(this);
		}
	}
	
	public List<CloneItem> getCloneInstances() {
		List<CloneItem> copyToReturn = new ArrayList<CloneItem>();
		for (CloneItem cloneItem : CloneItems)
			copyToReturn.add(cloneItem);
		return copyToReturn;
	}
	
	public void setTestPackages(String[] testPackages) {
		this.testPackages.clear();
		if (testPackages != null)
			for (String packageName : testPackages)
				this.testPackages.add(packageName);
	}
	
	public Set<String> getTestPackages() {
		return this.testPackages;
	}
	
	public void setTestSourceFolders(String[] testSourceFolders) {
		this.testSourceFolders.clear();
		if (testSourceFolders != null)
			for (String testSourceFolder : testSourceFolders)
				this.testSourceFolders.add(testSourceFolder);
	}
	
	public Set<String> getTestSourceFolders() {
		return this.testSourceFolders;
	}

}
