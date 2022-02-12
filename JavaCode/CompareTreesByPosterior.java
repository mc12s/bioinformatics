//This program trims the tree down to just topology (e.g. ((A,(B,C)),O); ) and determines if it matches the guide species tree
//usage: java IdentifyTrees T_10000_5000/1/r0/1024_1_25/u_.015 2 1M annotated_species.trees
import java.io.*;

public class CompareTreesByPosterior{

	public static void main(String[] args){

		int beast_reps = Integer.parseInt(args[1]);
		System.out.println(args[0]+"\t");
		String chainlength = args[2];
		String file_end = args[3];
		//boolean verbose = false;
		//if(args[3].equals("verbose")){
		//	verbose=true;
		//}

		//String[] filenames = {"5000_1.nwk", "5000_2.nwk"};

		int matches=0;

		try{
		for(int seqgenrep=1; seqgenrep<=10; seqgenrep++){

			String translated_tree_1 = "";
			for(int rep=1; rep<=beast_reps; rep++){
				char[] taxa_labels = new char[5];
				//System.out.print(rep +"\t");
				String file_name = args[0] + "/rep" + seqgenrep + "/beast_" + chainlength  + "_msout_" + rep + "/" + file_end;
				//String file_name = filenames[rep-1];
				File inFile = new File(file_name);

				BufferedReader br = new BufferedReader(new FileReader(inFile));
				String line = "";
				for(int i=0; i<13; i++){		//skip header lines
					line = br.readLine();
				}
				for(int i=1; i<=4; i++){			
					line = br.readLine();
					//System.out.println(line.charAt(5) + " = " + line.charAt(7));
					//taxa_labels[Integer.parseInt(line.substring(5,6))] = line.charAt(7);
				}

				for(int i=0; i<2; i++){		//skip header lines
					line = br.readLine();
				}
				//System.out.println(line);
				boolean skipbracket = false;
				boolean skiplength = false;
				String trimmed_tree = "";

				boolean found_start = false;

				//get the posteriors from the tree
				String[] posterior_split = line.split("posterior=");
				double[] posteriors = new double[3];
				int p_index = 0;

				for(String part : posterior_split){
					if (part == posterior_split[0]) continue;
					//System.out.println(part.substring(0,part.indexOf("]")));
					posteriors[p_index] = Double.parseDouble(part.substring(0,part.indexOf("]")));
					p_index++;
				}

				//for(double posterior : posteriors){
				//	System.out.println(posterior);
				//}
				

				for(int i=0; i<line.length(); i++){
					char cur = line.charAt(i);
					if(!found_start){	//skip until open parenthesis
						if(cur=='('){
							found_start=true;
						}else{
							continue;
						}
					}
					if(skipbracket){			//strip everything between brackets
						if(cur == ']'){
							skipbracket=false;
						}
					}else if(skiplength){
						if(cur == ',' || cur == ')'){
							//System.out.print(cur);
							trimmed_tree += cur;
							skiplength = false;
						}
					}else{
						if(cur == '['){
							skipbracket=true;
						}else if(cur == ':'){
							skiplength=true;
						}else{
							trimmed_tree += cur;
							//System.out.print(line.charAt(i));			
						}

					}
				}

				//System.out.println(trimmed_tree);

				//trimmed_tree = trimmed_tree.replace('1', taxa_labels[1]);
				//trimmed_tree = trimmed_tree.replace('2', taxa_labels[2]);
				//trimmed_tree = trimmed_tree.replace('3', taxa_labels[3]);
				//trimmed_tree = trimmed_tree.replace('4', taxa_labels[4]);

				//System.out.println(trimmed_tree);



				/*if(verbose){
					if(!trimmed_tree.substring(1,2).equals("O") && !trimmed_tree.substring(11,12).equals("O")){
						//System.out.print(trimmed_tree.substring(1,2));
						System.out.print("Outgroup misplaced");
						System.out.println("\t"+ trimmed_tree);
					}else{		//Outgroup correctly placed
						if(!trimmed_tree.contains("B,C") && !trimmed_tree.contains("C,B")){
							System.out.print("No B,C clade");
							System.out.println("\t"+ trimmed_tree);
						}else{
							System.out.println("Correct Tree");
							//System.out.println("((A,(B,C)),O)");
						}
					}	

				}*/

				String translated_tree = "";

				
				//System.out.println(trimmed_tree);

				int comma_index = 0;
				int posterior_index = 0;
				for(int i=0; i<3; i++){	//find the commas and see if they're surrounded by 2 taxa to find inner clade
					comma_index = trimmed_tree.indexOf(",", comma_index+1);					
					//System.out.println(comma_index);
					if(Character.isDigit(trimmed_tree.charAt(comma_index-1)) && Character.isDigit(trimmed_tree.charAt(comma_index+1))){
						posterior_index = i;
						System.out.print(posteriors[i]+"\t");
						break;
					}
				}

				//System.out.println("----------------------"+trimmed_tree);

				switch(trimmed_tree){
					case "((1,(2,3)),4)" : translated_tree = "((A,(B,C)),O)";break;
					case "((1,(3,2)),4)" : translated_tree = "((A,(B,C)),O)";break;
					case "(((2,3),1),4)" : translated_tree = "((A,(B,C)),O)";break;
					case "(((3,2),1),4)" : translated_tree = "((A,(B,C)),O)";break;
					case "(4,(1,(2,3)))" : translated_tree = "((A,(B,C)),O)";break;
					case "(4,(1,(3,2)))" : translated_tree = "((A,(B,C)),O)";break;
					case "(4,((2,3),1))" : translated_tree = "((A,(B,C)),O)";break;
					case "(4,((3,2),1))" : translated_tree = "((A,(B,C)),O)";break;

                    case "((2,(1,3)),4)" : translated_tree = "((B,(A,C)),O)";break;
                    case "((2,(3,1)),4)" : translated_tree = "((B,(A,C)),O)";break;
                    case "(((1,3),2),4)" : translated_tree = "((B,(A,C)),O)";break;
                    case "(((3,1),2),4)" : translated_tree = "((B,(A,C)),O)";break;
                    case "(4,(2,(1,3)))" : translated_tree = "((B,(A,C)),O)";break;
                    case "(4,(2,(3,1)))" : translated_tree = "((B,(A,C)),O)";break;
                    case "(4,((1,3),2))" : translated_tree = "((B,(A,C)),O)";break;
                    case "(4,((3,1),2))" : translated_tree = "((B,(A,C)),O)";break;

                    case "((3,(2,1)),4)" : translated_tree = "((C,(B,A)),O)";break;
                    case "((3,(1,2)),4)" : translated_tree = "((C,(B,A)),O)";break;
                    case "(((2,1),3),4)" : translated_tree = "((C,(B,A)),O)";break;
                    case "(((1,2),3),4)" : translated_tree = "((C,(B,A)),O)";break;
                    case "(4,(3,(2,1)))" : translated_tree = "((C,(B,A)),O)";break;
                    case "(4,(3,(1,2)))" : translated_tree = "((C,(B,A)),O)";break;
                    case "(4,((2,1),3))" : translated_tree = "((C,(B,A)),O)";break;
                    case "(4,((1,2),3))" : translated_tree = "((C,(B,A)),O)";break;

                    case "((4,(2,3)),1)" : translated_tree = "((O,(B,C)),A)";break;
                    case "((4,(3,2)),1)" : translated_tree = "((O,(B,C)),A)";break;
                    case "(((2,3),4),1)" : translated_tree = "((O,(B,C)),A)";break;
                    case "(((3,2),4),1)" : translated_tree = "((O,(B,C)),A)";break;
                    case "(1,(4,(2,3)))" : translated_tree = "((O,(B,C)),A)";break;
                    case "(1,(4,(3,2)))" : translated_tree = "((O,(B,C)),A)";break;
                    case "(1,((2,3),4))" : translated_tree = "((O,(B,C)),A)";break;
                    case "(1,((3,2),4))" : translated_tree = "((O,(B,C)),A)";break;

                    case "((4,(1,3)),2)" : translated_tree = "((O,(C,A)),B)"; break;
                    case "((4,(3,1)),2)" : translated_tree = "((O,(C,A)),B)"; break;
                    case "(((1,3),4),2)" : translated_tree = "((O,(C,A)),B)"; break;
                    case "(((3,1),4),2)" : translated_tree = "((O,(C,A)),B)"; break;
                    case "(2,(4,(2,1)))" : translated_tree = "((O,(C,A)),B)"; break;
                    case "(2,(4,(3,1)))" : translated_tree = "((O,(C,A)),B)"; break;
                    case "(2,((1,3),4))" : translated_tree = "((O,(C,A)),B)"; break;
                    case "(2,((3,1),4))" : translated_tree = "((O,(C,A)),B)"; break;

                    case "((4,(2,1)),3)" : translated_tree = "((O,(B,A)),C)"; break;
                    case "((4,(1,2)),3)" : translated_tree = "((O,(B,A)),C)"; break;
                    case "(((2,1),4),3)" : translated_tree = "((O,(B,A)),C)"; break;
                    case "(((1,2),4),3)" : translated_tree = "((O,(B,A)),C)"; break;
                    case "(3,(4,(2,1)))" : translated_tree = "((O,(B,A)),C)"; break;
                    case "(3,(4,(1,2)))" : translated_tree = "((O,(B,A)),C)"; break;
                    case "(3,((2,1),4))" : translated_tree = "((O,(B,A)),C)"; break;
                    case "(3,((1,2),4))" : translated_tree = "((O,(B,A)),C)"; break;

                    case "((1,(2,4)),3)" : translated_tree = "((A,(B,O)),C)"; break;
                    case "((1,(4,2)),3)" : translated_tree = "((A,(B,O)),C)"; break;
                    case "(((2,4),1),3)" : translated_tree = "((A,(B,O)),C)"; break;
                    case "(((4,2),1),3)" : translated_tree = "((A,(B,O)),C)"; break;
                    case "(3,(1,(2,4)))" : translated_tree = "((A,(B,O)),C)"; break;
                    case "(3,(1,(4,2)))" : translated_tree = "((A,(B,O)),C)"; break;
                    case "(3,((2,4),1))" : translated_tree = "((A,(B,O)),C)"; break;
                    case "(3,((4,2),1))" : translated_tree = "((A,(B,O)),C)"; break;

                    case "((3,(2,4)),1)" : translated_tree = "((C,(B,O)),A)"; break;
                    case "((3,(4,2)),1)" : translated_tree = "((C,(B,O)),A)"; break;
                    case "(((2,4),3),1)" : translated_tree = "((C,(B,O)),A)"; break;
                    case "(((4,2),3),1)" : translated_tree = "((C,(B,O)),A)"; break;
                    case "(1,(3,(2,4)))" : translated_tree = "((C,(B,O)),A)"; break;
                    case "(1,(3,(4,2)))" : translated_tree = "((C,(B,O)),A)"; break;
                    case "(1,((2,4),3))" : translated_tree = "((C,(B,O)),A)"; break;
                    case "(1,((4,2),3))" : translated_tree = "((C,(B,O)),A)"; break;

                    case "((2,(1,4)),3)" : translated_tree = "((B,(A,O)),C)"; break;
                    case "((2,(4,1)),3)" : translated_tree = "((B,(A,O)),C)"; break;
                    case "(((1,4),2),3)" : translated_tree = "((B,(A,O)),C)"; break;
                    case "(((4,1),2),3)" : translated_tree = "((B,(A,O)),C)"; break;
                    case "(3,(2,(1,4)))" : translated_tree = "((B,(A,O)),C)"; break;
                    case "(3,(2,(4,1)))" : translated_tree = "((B,(A,O)),C)"; break;
                    case "(3,((1,4),2))" : translated_tree = "((B,(A,O)),C)"; break;
                    case "(3,((4,1),2))" : translated_tree = "((B,(A,O)),C)"; break;

                    case "((3,(1,4)),2)" : translated_tree = "((C,(A,O)),B)"; break;
                    case "((3,(4,1)),2)" : translated_tree = "((C,(A,O)),B)"; break;
                    case "(((1,4),3),2)" : translated_tree = "((C,(A,O)),B)"; break;
                    case "(((4,1),3),2)" : translated_tree = "((C,(A,O)),B)"; break;
                    case "(2,(3,(1,4)))" : translated_tree = "((C,(A,O)),B)"; break;
                    case "(2,(3,(4,1)))" : translated_tree = "((C,(A,O)),B)"; break;
                    case "(2,((1,4),3))" : translated_tree = "((C,(A,O)),B)"; break;
                    case "(2,((4,1),3))" : translated_tree = "((C,(A,O)),B)"; break;

                    case "((1,(4,3)),2)" : translated_tree = "((A,(C,O)),B)"; break;
                    case "((1,(3,4)),2)" : translated_tree = "((A,(C,O)),B)"; break;
                    case "(((4,3),1),2)" : translated_tree = "((A,(C,O)),B)"; break;
                    case "(((3,4),1),2)" : translated_tree = "((A,(C,O)),B)"; break;
                    case "(2,(1,(4,3)))" : translated_tree = "((A,(C,O)),B)"; break;
                    case "(2,(1,(3,4)))" : translated_tree = "((A,(C,O)),B)"; break;
                    case "(2,((4,3),1))" : translated_tree = "((A,(C,O)),B)"; break;
                    case "(2,((3,4),1))" : translated_tree = "((A,(C,O)),B)"; break;

                    case "((2,(4,3)),1)" : translated_tree = "((B,(C,O)),A)"; break;
                    case "((2,(3,4)),1)" : translated_tree = "((B,(C,O)),A)"; break;
                    case "(((4,3),2),1)" : translated_tree = "((B,(C,O)),A)"; break;
                    case "(((3,4),2),1)" : translated_tree = "((B,(C,O)),A)"; break;
                    case "(1,(2,(4,3)))" : translated_tree = "((B,(C,O)),A)"; break;
                    case "(1,(2,(3,4)))" : translated_tree = "((B,(C,O)),A)"; break;
                    case "(1,((4,3),2))" : translated_tree = "((B,(C,O)),A)"; break;
                    case "(1,((3,4),2))" : translated_tree = "((B,(C,O)),A)"; break;

                    case "((1,4),(2,3))" : translated_tree = "((A,O),(B,C))"; break;
                    case "((1,4),(3,2))" : translated_tree = "((A,O),(B,C))"; break;
                    case "((4,1),(2,3))" : translated_tree = "((A,O),(B,C))"; break;
                    case "((4,1),(3,2))" : translated_tree = "((A,O),(B,C))"; break;
                    case "((2,3),(1,4))" : translated_tree = "((A,O),(B,C))"; break;
                    case "((2,3),(4,1))" : translated_tree = "((A,O),(B,C))"; break;
                    case "((3,2),(1,4))" : translated_tree = "((A,O),(B,C))"; break;
                    case "((3,2),(4,1))" : translated_tree = "((A,O),(B,C))"; break;

                    case "((2,4),(1,3))" : translated_tree = "((B,O),(A,C))"; break;
                    case "((2,4),(3,1))" : translated_tree = "((B,O),(A,C))"; break;
                    case "((4,2),(1,3))" : translated_tree = "((B,O),(A,C))"; break;
                    case "((4,2),(3,1))" : translated_tree = "((B,O),(A,C))"; break;
                    case "((1,3),(2,4))" : translated_tree = "((B,O),(A,C))"; break;
                    case "((1,3),(4,2))" : translated_tree = "((B,O),(A,C))"; break;
                    case "((3,1),(2,4))" : translated_tree = "((B,O),(A,C))"; break;
                    case "((3,1),(4,2))" : translated_tree = "((B,O),(A,C))"; break;

                    case "((3,4),(2,1))" : translated_tree = "((C,O),(B,A))"; break;
                    case "((3,4),(1,2))" : translated_tree = "((C,O),(B,A))"; break;
                    case "((4,3),(2,1))" : translated_tree = "((C,O),(B,A))"; break;
                    case "((4,3),(1,2))" : translated_tree = "((C,O),(B,A))"; break;
                    case "((1,2),(3,4))" : translated_tree = "((C,O),(B,A))"; break;
                    case "((1,2),(4,3))" : translated_tree = "((C,O),(B,A))"; break;
                    case "((2,1),(3,4))" : translated_tree = "((C,O),(B,A))"; break;
                    case "((2,1),(4,3))" : translated_tree = "((C,O),(B,A))"; break;

					default: System.out.println("No matching case found for " + trimmed_tree); break;
				}

                if(rep==1){
                	translated_tree_1 = translated_tree;
                }else{
                	if(translated_tree_1.equals(translated_tree)){
                	//System.out.println("Beast runs match");
					System.out.println("1");
					matches++;
                    }else{
                    //System.out.println(translated_tree_1+" mismatch "+translated_tree);
					System.out.println("0");
                	}
				}
				
				
				//System.out.println(translated_tree);
			}
		}
		}catch(IOException ioe){
			System.out.println(ioe);
		}
		System.out.println(matches);

	}

}
