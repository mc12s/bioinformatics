//BUG EXISTS species are not always given the right branch length totals depending
// on the treestring. Works fine for 3 species


import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class GuideTree{

	public static int[][] getSpeciationGroups(int nInds, String folder){
		
		String treestring = "";
		//System.out.println("________________________________________________"+nInds+"-------------------------");
		try{
			treestring = readTreeFile(folder + "guidetree.tre");
		}
		catch(IOException ioe){
			System.out.println("ERROR: Could not find Guide Tree file in " + folder +"guidetree.tre");
			System.exit(0);
		}

		//treestring = "((A:1000,B:1000):2500,(C:2000,(D:1500,E:1500,F:1500):500):1500):6500;";
		//treestring = "((A:1000,B:1000):2500,((D:1500,E:1500,F:1500):500,C:2000):1500):6500;";
		//treestring = "((B:1500,C:1500,D:1500):500,A:2000):1500;";
		//treestring = "(A:2000,(B:1500,C:1500,D:1500):500):1500;";
	
		int speciescount=1;
		for(int i=0; i<treestring.length(); i++){
			if(treestring.charAt(i)==','){
				speciescount++;
			}
		}

		Tree T = new Tree(0, null, treestring, 1, nInds);
		
		getClade(treestring, T, nInds);
		//System.out.println("-------------------");
		//System.out.println();

		ArrayList<String> speciationlist = new ArrayList<String>();
		T.traverse(T, 0, speciationlist);

		//System.out.println(speciationlist + "\n");

		//System.out.println("Unsorted list");

		//for(String entry : speciationlist){
		//	System.out.println("Entry: " + entry);
		//}


		String[] speciationArray = new String[speciationlist.size()];

		Collections.sort(speciationlist, new Comparator<String>(){	//overriding sort to compare speciation list entries by their first string as an integer
			public int compare(String s1, String s2){
					Integer val1 = Integer.parseInt(s1.split(" ")[0]);
					Integer val2 = Integer.parseInt(s2.split(" ")[0]);
					return val1.compareTo(val2);
				}
		});

		//Collections.sort(speciationlist, string_comparator);
		int entryindex = 0;

		//System.out.println("Sorted list");


		//for(String entry : speciationlist){
		//	System.out.println("Entry: " + entry);
		//}


		for(String entry : speciationlist){
			if(entry.charAt(0) != '0'){
				String[] parts = entry.split(" ");
				String newentry = parts[0];
					if(Integer.parseInt(parts[1])>1){
						for(int i=0;i<entryindex;i++){
							newentry+=","+String.valueOf(nInds);
						}
						for(int i=2; i<parts.length;i++){
							newentry+="," + parts[i];
						}
					}
					else{
						for(int i=2; i<parts.length;i++){
							newentry+="," + parts[i];
						}
						for(int i=0;i<entryindex;i++){
							newentry+=","+String.valueOf(nInds);
						}				
					}

				speciationArray[entryindex] = newentry;
				entryindex++;
			}
		}

		int colsize = 2*(speciationArray[entryindex-1].split(",").length)-1;
		int rowsize = speciationlist.size();
		int speciationGroups[][] = new int[rowsize*10][colsize];

		int row_iter=0;
		for(String entry : speciationArray){
			int cur_total = 0;
			int increment = 0;
			//System.out.println("\nentry: " + entry);
			String parts[];
			if(entry!=null){
				parts = entry.split(",");
			}
			else{
				break;
			}
			int part_iter = 0;
			int inds_total = 0;
			for(String part : parts){
				if(part_iter==0){
					speciationGroups[row_iter][part_iter] = Integer.parseInt(part);
				}			
				else{
					//System.out.print(Integer.parseInt(part) + ": ");


					speciationGroups[row_iter][increment+part_iter] = (Integer.parseInt(part)-Integer.parseInt(part)+cur_total);
					speciationGroups[row_iter][increment+part_iter+1] = (Integer.parseInt(part)+cur_total-1);

					increment+=1;

					//System.out.print((Integer.parseInt(part)-Integer.parseInt(part)+cur_total) + " " + (Integer.parseInt(part)+cur_total-1) + " ");
					cur_total+=Integer.parseInt(part);
				}
				part_iter++; 
			}
			row_iter++;
		}

		//System.out.println();
		speciationGroups[entryindex][0] = -999999;

		//for(String entry : speciationArray){
		//	System.out.println(entry);
		//}

		//System.out.println();

		//System.out.println("Checkpoint");

		//for(int i=0; i<rowsize; i++){
		//	for(int j=0; j<colsize; j++){
		//		System.out.print(speciationGroups[i][j] + " ");
		//	}
		//	System.out.println();
		//}

		return speciationGroups;
	}




	public static String getClade(String cladestring, Tree T, int nInds){
		int openparen=0;
		int firstopen=-1;
		int species=1;
		int speciespos=1;
		String branchstring = "";

		for(int i=0; i<cladestring.length(); i++){
			if(cladestring.charAt(i)=='('){
				openparen++;
				if(firstopen==-1){
					firstopen=i;
				}
			}
			else if(cladestring.charAt(i)==')'){
				openparen--;
				if(openparen==0){
					//System.out.println(i);
					branchstring = "1";

					int branchlength = getBranchlength(cladestring.substring(i+1));

					String newcladestring=cladestring.substring(firstopen+1,i);
					System.out.println(newcladestring + " Species: " + species);		

					System.out.println("parent branchlength = " + T.branchlength + " + " + branchlength);
					Tree node = new Tree(branchlength + T.branchlength, T, newcladestring, speciespos, nInds);
					T.leaves.add(node);

					getClade(newcladestring, node, nInds);
					firstopen=-1;
				}
			}
			else if(firstopen==-1 && cladestring.charAt(i)==','){
				species++;
				speciespos=species;
				//String[] parts = cladestring.split(",");
				//for(String part : parts){
				//	System.out.print(part + " ");
				//}
				//System.out.println();
			}
			else if(i==cladestring.length()-1){
				T.setSpecies(species);
			}
		}
		return cladestring;

	}

	public static int getBranchlength(String selection){
		//System.out.println("getBranchlength : " + selection);
		String branchstring = "";
		for(int i=1; i<selection.length(); i++){
			if(Character.isDigit(selection.charAt(i))){
				branchstring+=selection.charAt(i);
			}
			else{
				break;
			}
		}
		int branchlength = Integer.parseInt(branchstring);
		return branchlength;
	}


	public static int[] splitPop(int nInds, int divisions){
	//This method splits an int into d parts evenly then add the remainder across the parts
	// e.g. 20/3 = 7 7 6
		//System.out.println("inside splitPop\n\tspecies=" + divisions + "\n");
		if(divisions==0){ int parts[] = {nInds}; return parts;}
		
		int n = nInds;
		int d = divisions;
		int[] parts = new int[d];
		//System.out.println(d);

		int m=n/d;
		for(int i=0; i<d; i++){
			parts[i] = m;		
		}

		int r=n%d;
		for(int i=0;i<r;i++){
			parts[i] += 1;
		}

		return parts;		

		//for(int i=0; i<d; i++){
		//	System.out.print("\t"+parts[i]);		
		//}	
	}


	public static String readTreeFile(String filein) throws IOException{

		BufferedReader br = new BufferedReader(new FileReader(filein));

		String treestring = br.readLine();

		br.close();

		return treestring;

	}

}

class Tree{

	//String treestring = "((A:1000,B:1000):2500,(C:2000,(D:1500,E:1500,F:1500):500):1500):6500;";
	
	//int[][] speciationGroups = {{5000, 0,9, 10,19},{7300, 0,9, 10,19, 20,39},{8000, 0,19, 20,39, 40,49, 50,59},{9000, 0,19, 20,39, 40,59, 60,69, 70,79},{-999999}};

	//int[][] speciationGroups = {{6500, 0,9, 10,19},{8000, 0,19, 20, 29, 30, 39},{8500, 0,19, 20,39, 40,46, 47,53, 54,59},{9000, 0,9, 10,19, 20,39, 40,59, 60,79, 80,99},{-999999}};

	Tree parent = null;
	int branchlength; 
	ArrayList<Tree> leaves;
	private int nInds;
	String nodestring="";
	private int species;
	private int speciespos;


	public Tree(int branchlength, Tree parent, String nodestring, int speciespos, int nInds){
		this.nInds = nInds;
		this.branchlength = branchlength;
		this.parent = parent;
		this.nodestring = nodestring;
		this.speciespos = speciespos;
		this.leaves = new ArrayList<Tree>();
		//System.out.println("constructor: " + branchlength);
		//System.out.println(species);
		//System.out.println();
	}

	public void setSpecies(int species){
		this.species=species;
	}

	public void traverse(Tree node, int depth, ArrayList<String> speciationlist){
		//System.out.println("generation: " + node.branchlength);
		//System.out.println("clade: " + node.nodestring);
		//System.out.println("Species: " + node.species);
		//System.out.println("depth = " + depth + " leaves = " + node.leaves.size());
		depth++;

		int[] parts = GuideTree.splitPop(nInds, node.species );
		String speciationString = node.branchlength + " " + node.speciespos +  " ";

		for(int part : parts){
			speciationString += (part + " ");
		}

		speciationlist.add(speciationString);
		//System.out.println("\n");
		for(int i=0; i<node.leaves.size();i++){
			//System.out.println("depth = " + depth + " leaf: " + (i+1));

			if(node.leaves.get(i) != null){
				traverse (node.leaves.get(i), depth, speciationlist);
			}
			else{
				depth--;
			}
		}
	}
}
