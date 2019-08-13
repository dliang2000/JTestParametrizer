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
	
	private long subtreeMatchingTimeSum = 0;
	private int numberOfStatementsToBeRefactoredSum = 0;
	private int numberOfNodeComparisonsSum = 0;
	private long subtreeMatchingWallNanoTimeSum = 0;
	private AnalysisStatus status;
	private List<List<PDGSubTreeMapperInfo>> subtreeMappersListList = new ArrayList<List<PDGSubTreeMapperInfo>>();
	private Set<String> testPackages = new HashSet<>();
	private Set<String> testSourceFolders = new HashSet<>();

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
	
	public void setSubtreeMatchingTimeSum(long[] subtreeMatchingTimeList) {
		this.subtreeMatchingTimeSum = 0;
		if (subtreeMatchingTimeList != null)
			for (long subtreeMatchingTime : subtreeMatchingTimeList)
				this.subtreeMatchingTimeSum +=  subtreeMatchingTime;
	}
	
	public long getSubtreeMatchingTimeSum() {
		return this.subtreeMatchingTimeSum;
	}
	
	public void setNumberOfStatementsToBeRefactoredSum(int[] numberOfStatementsToBeRefactoredList) {
		this.numberOfStatementsToBeRefactoredSum = 0;
		if (numberOfStatementsToBeRefactoredList != null)
			for (int numberOfStatementsToBeRefactored : numberOfStatementsToBeRefactoredList)
				this.numberOfStatementsToBeRefactoredSum += numberOfStatementsToBeRefactored;
	}
	
	public int getNumberOfStatementsToBeRefactoredSum() {
		return this.numberOfStatementsToBeRefactoredSum;
	}
	
	public void setNumberOfNodeComparisonsSum(int[] numberOfNodeComparisonsList) {
		this.numberOfNodeComparisonsSum = 0;
		if (numberOfNodeComparisonsList != null)
			for (int numberOfNodeComparisons : numberOfNodeComparisonsList)
				this.numberOfNodeComparisonsSum += numberOfNodeComparisons;
	}
	
	public int getNumberOfNodeComparisonsSum() {
		return this.numberOfNodeComparisonsSum;
	}
	
	public void setSubtreeMatchingWallNanoTimeSum(long[] subtreeMatchingWallNanoTimeList) {
		this.subtreeMatchingWallNanoTimeSum = 0;
		if (subtreeMatchingWallNanoTimeList != null)
			for (long subtreeMatchingWallNanoTime : subtreeMatchingWallNanoTimeList)
				this.subtreeMatchingWallNanoTimeSum += subtreeMatchingWallNanoTime;
	}
	
	public long getSubtreeMatchingWallNanoTimeSum() {
		return this.subtreeMatchingWallNanoTimeSum;
	}
	
	public AnalysisStatus getStatus() {
		return status;
	}

	public void setStatus(AnalysisStatus status) {
		this.status = status;
	}
	
	public List<List<PDGSubTreeMapperInfo>> getPDGSubTreeMappersInfoListList() {
		return this.subtreeMappersListList;
	}
	
	public void clearMappersInfo() {
		this.subtreeMappersListList = new ArrayList<List<PDGSubTreeMapperInfo>>();
	}
	
	public boolean getRefactorable() {
		int count = 0;
	  	for (List<PDGSubTreeMapperInfo> pdgSubTreeList: subtreeMappersListList) {
	  		for (PDGSubTreeMapperInfo info: pdgSubTreeList) {
	  			if(info.getRefactoringWasOK()) {
	  				count++;
	  				break;
	  			}
	  		}
	 	}
	  	return count == cloneGroupSize - 1;
	}
	
	/**
	 * Get the list of only refactorable mappers
	 * @return
	 */
	public List<List<PDGSubTreeMapperInfo>> getRefactorableMappersInfo() {
		List<List<PDGSubTreeMapperInfo>> toReturn = new ArrayList<List<PDGSubTreeMapperInfo>>();
		for (List<PDGSubTreeMapperInfo> pdgSubTreeList : subtreeMappersListList) {
			List<PDGSubTreeMapperInfo> pdgSubTreeToAdd = new ArrayList<>();
			for (PDGSubTreeMapperInfo info: pdgSubTreeList) {
	  			if(info.getRefactoringWasOK()) {
	  				pdgSubTreeToAdd.add(info);
	  			}
	  		}
			toReturn.add(pdgSubTreeToAdd);
		}
		return toReturn;
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
