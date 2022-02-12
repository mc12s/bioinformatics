import java.util.LinkedList;

public class Tree{

	Tree parent = null;
	private int branchlength; 
	private LinkedList<Tree> leaves;
	private int nInds = 20;

	public Tree(int branchlength, Tree parent){
		this.branchlength = branchlength;
		this.parent = parent;
		leaves = new LinkedList<Tree>();
		System.out.println("constructor: " + branchlength);
	}

	public void traverse(Tree node, int depth){
		System.out.println("branch length: " + node.branchlength);
		System.out.println("depth = " + depth + " leaves = " + node.leaves.size());
		depth++;
		int[] parts = splitPop(nInds, node.leaves.size()+1 );
		for(int part : parts){
			System.out.print(part + " ");
		}
		System.out.println("\n");
		for(int i=0; i<node.leaves.size();i++){
			System.out.println("depth = " + depth + " leaf: " + (i+1));
			if(node.leaves.get(i) != null){
				traverse (node.leaves.get(i), depth);
			}
			else{
				depth--;
			}
		}

	}

	public static void main(String args[]){
		
		int nInds = 20;

		String treestring = "((A:1000,B:1000):2500,(C:2000,(D:1500,E:1500):500):1500):6500;";
		int species=1;

		int branchlength = Integer.parseInt( treestring.substring( treestring.lastIndexOf(":")+1, treestring.length()-1 ));

		Tree T = new Tree(branchlength, null);
		
		for(int i=0; i<treestring.length()-1; i++){
			if(treestring.charAt(i) == ','){
				species++;
			}
		}
		getClade(treestring, T);
		System.out.println();
		T.traverse(T, 0);
	}


	public static String getClade(String cladestring, Tree T){
		int first_open = 0;
		int open_paren = 0;

		String newcladestring = "";
		String branchstring = "";

		for(int i=1; i<cladestring.length()-1; i++){
			if(cladestring.charAt(i) == '('){
				open_paren++;
				if(open_paren==1){first_open=i;};				
			}
			else if(cladestring.charAt(i) == ')'){
				open_paren--;
				if(open_paren==0){ 
					branchstring = "";
					for(int j=2; j<cladestring.length(); j++){
						int k=i+j;
						if(Character.isDigit(cladestring.charAt(k))){
							branchstring+=cladestring.charAt(k);
						}
						else{
							break;
						}
					}
					int branchlength = Integer.parseInt(branchstring);
					System.out.println(i + " " + cladestring + " " + Integer.parseInt(branchstring));				
					newcladestring = cladestring.substring(first_open,i+1);

					System.out.println("parent branchlength = " + T.branchlength);
					Tree node = new Tree(branchlength+T.branchlength,T);
					T.leaves.add(node);
					getClade(newcladestring, node);
				}
			}
			branchstring="";
		}
		return newcladestring;
	}

	public static int[] splitPop(int nInds, int divisions){
	//This method splits an int into d parts evenly then add the remainder across the parts
	// e.g. 20/3 = 7 7 6

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

}
