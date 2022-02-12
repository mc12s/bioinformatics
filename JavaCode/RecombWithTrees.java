//Usage: java RecombWithTrees nInds gens chunks recomb_rate speciation_time taxasamples threads foldername replicate
//e.g. java RecombWithTrees 10 100 5 1 100
//java -Xms4g -Xmx13g RecombWithTrees 10000 500000 2 0 8 8 T_5000_2500/ 1 false
//java -Xmx10g -Xms10g RecombWithTrees 10000 97500 2 .5 8 8 T_5000_2500/ 9999 true

import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;
import java.util.Random;
import java.util.HashMap;
import java.util.Set;
import java.util.Iterator;
import java.util.Map;
import java.util.Date;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;

public class RecombWithTrees{

	//public static int isCommaCounter=0;
	//public static int isNoLengthCounter=0;

	public static void main(String args[]){
	boolean loop_until_crash = false;
	do{ //loop program until it crashes

		Date date = new Date();
		System.out.println(date.toString());

		String treestring = "";
		int nInds=Integer.parseInt(args[0]);
		int gens=Integer.parseInt(args[1]);
		int numNuc=Integer.parseInt(args[2]);
		double recombRateAll=Double.parseDouble(args[3]);
		double[] recombRate = new double[numNuc];
		int selected_species = 3;
		int nTips=Integer.parseInt(args[4]);
		String foldername = args[6];
		String replicate = args[7];
		//boolean burn_in_check = Boolean.parseBoolean(args[8]);
		boolean burn_in_check = false;
		if(burn_in_check){ nTips = 33; }



		Random r1 = new Random();

		/*int[][] selectedTips = {
				{(int)(.1*nInds)},
				{(int)(1.1*nInds)},
				{(int)(2.1*nInds)}
		};
		*/	

		int[][] selectedTips = chooseRandomTips(nInds,selected_species,nTips);
	
		for(int i=0; i<selected_species; i++){
			for(int entry:selectedTips[i]){
				System.out.print(entry + " ");
			}
			System.out.println();
		}

	  ArrayList<Chunk> SelectedChunkArr = new ArrayList<Chunk>();		
		ArrayList<Ind> CurGenArr = new ArrayList<Ind>(nInds);	

		final int NUMTHREADS=Integer.parseInt(args[5]);

		ThreadPoolExecutor executor = (ThreadPoolExecutor) Executors.newFixedThreadPool(NUMTHREADS);
		//ThreadPoolExecutor executor = (ThreadPoolExecutor) Executors.newCachedThreadPool();

		//if(NUMTHREADS<1){NUMTHREADS=1;}

		for(int i=0; i<recombRate.length;i++){													//all recombination rates for different chunks are the same for now
			recombRate[i]=recombRateAll;
		}

		//Chunk rootChunk = new Chunk(null);

		//for(int i=0; i<nInds; i++) {
		//	RootGenArr.add(new Ind(rootChunk) );
		//}

		for(int i=0; i<nInds; i++){
			CurGenArr.add(new Ind(numNuc) );
		}

		//for(int i=0; i<nInds; i++){
		//	CurGenArr.add(RootGenArr.get(i));
		//}

		//int[] ancestors_id = new int[2];
		
		//int NUMTHREADS = 8;

		//new generation

		int species = 1;																									//species is how many species exist (speciation events + 1)
		int[][] speciationArray = getSpeciationArray(nInds, foldername);


		System.out.println("SpeciationArray:");
		
		for(int j=0;j<3;j++){	
			for(int i=0;i<speciationArray[j].length;i++){
				System.out.print("\t" + speciationArray[j][i] + " ");
			}		
			System.out.println();		
		}
		System.out.println();		
		
		int[] speciationSubArray= new int[10];
		speciationSubArray=speciationArray[0];
		int speciation_time=speciationSubArray[0];
		boolean speciation = false;


		int[] splits = splitInds(nInds, NUMTHREADS);
		int[] splitadd = new int[splits.length];

		for(int i=0; i<splitadd.length; i++){
			if(splits[i]%2==0){
				splitadd[i]=1;
			}
			else{
				splitadd[i]=2;
			}
		}



		for(int curGen=0;curGen<gens;curGen++){
			//int[] randomArray = new Random().ints(2*nInds, 0, nInds).toArray();
			//if(curGen%1000==0){	
				//System.out.print("gen: " + curGen + "\n");// + "\r");
				//System.out.println("Mem(while generating inds): " + ((Runtime.getRuntime().totalMemory() - 	Runtime.getRuntime().freeMemory())/1024/1024)+"Mb" );
				//date = new Date();
				//System.out.println(date.toString());
			//}


			//if(curGen==speciation_time){  System.out.println("speciation event:"+speciation_time); speciationArray = getSpeciationArray(species-1, nInds); speciation_time=getSpeciationArray(species, nInds)[0][0]; species++; speciation=true;}
			if(curGen==speciation_time){  System.out.print("\n\n---------------------\nspeciation event:"+speciation_time +"\n---------------------\n"); speciationSubArray = speciationArray[species-1]; speciation_time=speciationArray[species][0]; species++; speciation=true;}
			ArrayList<Ind> NewGenArr = new ArrayList<Ind>(species*nInds);	

			for(int j=0; j<species; j++){																	//for each new species, create nInds doubling but splitting the population for new gen
				int min=j*nInds;																						//normally nInds ancestors will have nInds children
				int max=j*nInds+nInds-1;
				if(speciation){
					min=speciationSubArray[j*2+1];																//... but during speciation we get the start and end ancestors from the speciationArray
					max=speciationSubArray[j*2+2];
					System.out.println("min="+min+" , max="+max);
				}
				if((curGen%1000)==1){System.out.print("\rgen:"+curGen + " min="+min+" , max="+max);}


				/*for(int i=0; i<nInds; i++) {
					ancestors_id[0]=randomMating(min,max,-1);
					ancestors_id[1]=randomMating(min,max,ancestors_id[0]);

					//System.out.print("curGen: " + curGen + " , species: " + (j+1) + " , parents: " +ancestors_id[0] + " " + ancestors_id[1] + "\r");
					
					NewGenArr.add(new Ind(ancestors_id, numNuc, recombRate) );
				}*/



 				//long start = System.nanoTime();
				int splitcount=0;

				ArrayList<NewThread> callables = new ArrayList<NewThread>(NUMTHREADS);
//			Thread[] runThreads = new Thread[NUMTHREADS];

				for(int i=1;i<=NUMTHREADS;i++){
					splitcount+=splits[i-1];
					callables.add(new NewThread(i, splitcount, splitcount+splits[i]-1, splitadd[i], min, max, numNuc, recombRate, CurGenArr));
				}

				try{ 
					List<Future<ArrayList<Ind>>> resultList = executor.invokeAll(callables);
						try{
				      for(Future<ArrayList<Ind>> result : resultList)
				 	     {
								NewGenArr.addAll(result.get());
				 	     }
							} catch (ExecutionException exex){
								exex.printStackTrace();
							} 
		        } catch (InterruptedException inex){
						inex.printStackTrace();
					}         

/*					splitcount+=splits[i-1];
					ThreadList[i-1]=new NewThread(i, splitcount, splitcount+splits[i]-1, splitadd[i], min, max, numNuc, recombRate, CurGenArr);
					//System.out.println("Splits " + splitcount + ", " + (splitcount+splits[i]-1));
					runThreads[i-1]=new Thread(ThreadList[i-1], "Thread_"+i);
					runThreads[i-1].start();
				}
*/

/*				for(int i=1;i<=runThreads.length;i++){
					try{
						//ThreadList.get(i-1).t.join();
						runThreads[i-1].join();
						NewGenArr.addAll(ThreadList[i-1].getThreadArray());
						//System.out.println(i);
					}catch(InterruptedException e){
						System.out.println("join interrupted");
					}

				}
*/
				//long time = System.nanoTime() - start;
		    //System.out.println("time all " + (double)time + "\n");


 				//start = System.nanoTime();
				//for(int i=1;i<=runThreads.length;i++){
				//		NewGenArr.addAll(ThreadList[i-1].getThreadArray());
				//}
				//time = System.nanoTime() - start;
		    //System.out.println("time merge " + (double)time);
				//System.out.println();
								
				//System.out.println("joined");
			}
		
			CurGenArr=NewGenArr;
			//CurGenArr.clear();
			//for(int i=0; i<nInds; i++) {
			//CurGenArr.addAll(NewGenArr);
			//}
			//NewGenArr.clear();
			NewGenArr=null;
			speciation=false;
		}

    executor.shutdown(); //kills the thread pool service
	
		//System.out.println("Mem(after generating inds): " + ((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024)+"Mb" );

		date = new Date();
		System.out.print("\n" + date.toString() + "\n");


		String allele_set = "0 1 r b";
		String[] allele_parts = allele_set.split(" ");
		ArrayList<Integer> SelectedTipNames = new ArrayList<Integer>();	

		int subsample_termination=1;
		if(burn_in_check){ //if checking burn in times select a random 100 tips and create a tree
			selectedTips = chooseRandomTips(nInds,selected_species,nTips);
			subsample_termination=selectedTips[0].length;
		}

		for(int burn_in_run=0; burn_in_run<subsample_termination; burn_in_run++){
			for(int subsample=selectedTips[0].length; subsample>=subsample_termination; subsample--){

				for(String allele_part : allele_parts){

				int[][] random_allele_choice = new int[selected_species][nTips];
				int[] random_allele_choice1d = new int[(selected_species*nTips)];
					if(allele_part.equals("r")){
						//System.out.println("YES");
						Random r_allele = new Random();
						for(int m=0; m<selected_species; m++){
							for(int i=0; i<nTips; i++){
								random_allele_choice[m][i] = r_allele.nextInt(2);
								//System.out.println(random_allele_choice[m][i]);
								random_allele_choice1d[(m*nTips)+i] = random_allele_choice[m][i];					//Store the random array in a 1d array so it's easier to use in the getTree method
							}
						}
					}

					for(int thisChunk=0; thisChunk<numNuc; thisChunk++){																															//we get a gene tree for each nucleotide chunk (exon)
							System.out.println("XXXXXXXXXXXXXXXXXXXXXXXXX " + CurGenArr.size() + " XXXXXXXXXXXXXXXXXXXXXXXXXX " + subsample + allele_part);
							//for(int i=0; i<CurGenArr.size(); i+=0){																									//this is the array of tips of the tree (the final generated generation) we start here and work backwards to track parentage
							for(int m=0; m<selected_species; m++){
								int random_iter=0;
								for(int i:selectedTips[m]){																																//we can specify the individuals we want to include in our tree here, selectedTips is user input
									String allele_choice = "";
									if(allele_part.equals("r")){						
										Random r_allele = new Random();
										//System.out.println("random_allele_choice array size: " + random_allele_choice[0].length() );
										allele_choice = Integer.toString(random_allele_choice[m][random_iter]);
										random_iter++;
										//System.out.println(allele_choice);
									}							
									else if(allele_part.equals("b")){
										allele_choice = "01";	
									}
									else{
										allele_choice = allele_part;
									}


									for(int allele=0; allele<allele_choice.length(); allele++){		//2 alleles
										if(i==-1){
											SelectedChunkArr.add(null);
											SelectedTipNames.add(null);
										}
										else{
											SelectedChunkArr.add(CurGenArr.get(i).Chrom[thisChunk][Integer.parseInt(allele_choice.substring(allele, allele+1))]);																	//this is the array that gets passed into getTree()
											SelectedTipNames.add(i);
										}
									}
								//i+=nInds/2;																																							//here we were skipping 10 inds each loop which results in 2 tip per species
								//																																												//POSSIBLE ERROR DUE TO DIVIDING an odd nInds and having a remainder
								}
							}


							HashMap<Chunk,String> pointer_map = new HashMap<Chunk,String>();
							//System.out.println("YYYYYYYYYYYYYYYYYYYYYYYYY " + SelectedChunkArr.size() + " YYYYYYYYYYYYYYYYYYYYYY ");
							boolean coalesce = getTree(pointer_map, SelectedChunkArr, SelectedTipNames, gens, allele_part, random_allele_choice1d);
				//UNCOMMENT FOR HASHMAP PRINTOUT
				//			Set set = pointer_map.entrySet();
				//			Iterator iter = set.iterator();

				//			while(iter.hasNext() ){
				//				Map.Entry entry = (Map.Entry)iter.next();
				//			System.out.print(entry.getKey() + ": ");
				//			System.out.println(entry.getValue() );
				//			}		
				//UNCOMMENT FOR TREE OUTPUT
							//System.out.println("Chunk " + j + ": " + pointer_map.get(null));
							//System.out.println("Chunk " + j + ": ");

							//treestring+=pointer_map.get(null)+";\n";
							if(coalesce){
								//treestring += "n" + String.valueOf(subsample) + "_c" + String.valueOf(thisChunk) + "_a" + allele_part + " (" + pointer_map.get(null) + ",n" + String.valueOf(subsample) + "_c" + String.valueOf(thisChunk) + "_a" + allele_part + ":1);\n";
								treestring += "n" + String.valueOf(subsample) + "_c" + String.valueOf(thisChunk) + "_a" + allele_part + " (" + pointer_map.get(null) + ",n" + String.valueOf(subsample) + "_cx_a" + allele_part + ":1);\n";
							}
							else{
								treestring += "n" + String.valueOf(subsample) + "_c" + String.valueOf(thisChunk) + "_a" + allele_part + "_Did_not_Coalesce\n"; 
							}

				//UNCOMMENT FOR MEMORY OUTPUT			
				//			System.out.println("Mem(for each tree): " + ((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024)+"Mb" );
							pointer_map=null;
							SelectedChunkArr.clear();
							SelectedTipNames.clear();
						}//end of chunk loop
						
					}//end of allele choice loop

					//This section loops through the selectedTips creating an array of 1 random less tip per row until we reach 1 selected tip.
					for(int m=0; m<selected_species; m++){
						//System.out.println("selectedtips[m] Length: " + selectedTips[m].length);
						int randomskip = randomTaxonSkip(selectedTips[m].length, selectedTips[m]);	//call to recursive method that sets a random index of selectedTips to -1
						//System.out.println("Randomskip=" + randomskip);
						selectedTips[m][randomskip] = -1;
					}

					for(int m=0; m<selected_species; m++){
						for(int ii=0; ii<selectedTips[m].length; ii++){
							System.out.print(" " + selectedTips[m][ii]);
						}
						System.out.println();
					}
					/*int[][] tempTips = new int[selected_species][subsample]; 
					Random r = new Random();
					for(int m=0; m<selected_species; m++){
						int randomskip = r.nextInt(subsample+1);
						int index=0;
						for(int i=0;i<subsample;i++){
							if(i!=randomskip){
								tempTips[m][index] = selectedTips[m][i];
								index++;
							}
						}
					}
					selectedTips=tempTips;
					*/

			}//end of subsampling

			if(burn_in_check){
				selectedTips = chooseRandomTips(nInds,selected_species,nTips); //get a new n tips
			}
		}

		treestring = treestring.replaceAll("\\)",".0\\)"); //regex to add .0 after branch length
		treestring = treestring.replaceAll(",",".0,"); 
		writeToFile(treestring, nInds, gens, numNuc, recombRateAll, speciationArray, replicate); 

		//System.out.println("no lengths=" + isNoLengthCounter);
		//System.out.println("comma between clades=" + isCommaCounter);
		date = new Date();
		System.out.println(date.toString());



	}while(loop_until_crash);
	}//end of main

	static int randomTaxonSkip(int subsample, int[] selectedTipsRow){  //recursive function that repeats a random int call if that index of selectedtips[m] is already -1
		Random r = new Random();
		int randomskip = r.nextInt(subsample);
		if(selectedTipsRow[randomskip]==-1){
				randomskip=randomTaxonSkip(subsample, selectedTipsRow);
		}
		return randomskip;		
	}

	static boolean getTree(HashMap<Chunk,String> pointer_map, ArrayList<Chunk>SelectedChunkArr, ArrayList<Integer>SelectedTipNames, int gens, String allele_part, int[] random_allele_choice1d){
		boolean coalesce = true;
		int arr_Iter=0;

		int allele_count = 1;
		if(allele_part.equals("b")){
			allele_count=2;
		}

		String new_tips = "";
		String tips = "AABBCCDDEEFFGGHHIIJJKKLLMMNNOOPPQQRRSSTTUUVVWWXXYYZZ";		//replace with taxon_set names in array
		for(int i=0; i<SelectedTipNames.size()/52;i++){
			new_tips +=tips;
		}
		tips+=new_tips;
		String[] taxon_set_names;

		//String cur_tip=tips.substring(1,2);
		String cur_tip=tips.substring(0,1);

		//System.out.println("SelectedChunkArray size: " + SelectedChunkArr.size());

		while(SelectedChunkArr.get(arr_Iter)==null){	//skip initial null subsample tips
		 arr_Iter+=allele_count;	
		}
		
		Chunk pChunk=SelectedChunkArr.get(arr_Iter).parent;
		

		String clade_val = "";
		String new_clade_val = "";
		String subname="0";

		int selectedArrSize = SelectedChunkArr.size();

		int mod1000 = 1000;

		while(arr_Iter<SelectedChunkArr.size() && pChunk!=null){				//checks tree 1 whole generation (rather than 1 ind's history) at a time
			if(SelectedChunkArr.get(arr_Iter) == null){
				arr_Iter+=allele_count;
				continue;																										//skip arr_Iter if this tip in the subsample was skipped
			}
			if(arr_Iter%mod1000==0){																
				int mapsize = pointer_map.size();
				boolean once = false;
				if(mapsize==1 && (SelectedChunkArr.get(arr_Iter).parent!=null)){																//if we've collapsed the tree to 1 string...
					if(mod1000==1000){		//this skips to the next generation to make sure clade_val gets wrapped in it's proper parenthesis and given a length
						mod1000=1;
					}
					else{					
					
						System.out.print("hashmap size " + mapsize +"        \r");
						Chunk sChunk = SelectedChunkArr.get(arr_Iter);								
						clade_val=pointer_map.get(sChunk);						//this is the whole string which needs the remaining length to gen 1 added.  
						//if(!once){	System.out.println(clade_val); once=true;}
						int start_length=clade_val.lastIndexOf(":")+1;
						int end_length=clade_val.substring(start_length).lastIndexOf(")"); //returns -1 if ")" does not exist
						//while(end_length<=clade_val.length() && Character.isDigit(clade_val.charAt(end_length)) ){ //This checks to find if there's a parenthesis at the end of the length, which I think is due to the double parenthesis wrap bug
						//	end_length++;
						//}
						int length=0;
						try{
							if(end_length!=-1){	//if clade_val ends in a ")" Theoretically should never arrive here anymore
								System.out.println("VVVVVVV-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------VVVVVVV");
								System.out.println(clade_val);
								System.out.println(clade_val.substring(start_length));
								System.out.println(clade_val.substring(start_length+end_length));
								System.exit(0);		
								length=Integer.parseInt(clade_val.substring(start_length,start_length+end_length));				//ERROR IS HERE, in a RARE case where two clades are joined on the final generation... fixed? yes
								if(!clade_val.substring(start_length+end_length).equals(")")){
									System.out.println("Error, exiting...");
									System.exit(0);							
								}
							}
							else{
								length=Integer.parseInt(clade_val.substring(start_length));
							}
						}catch(NumberFormatException | StringIndexOutOfBoundsException ex){
							System.out.println("Error:\t Start Length= "+start_length);
							System.out.println("\tEnd Length= "+end_length);
							System.out.println(ex.getMessage() );
							System.exit(0);
						}
						int newlength=0;
						while(sChunk!=null){																			//...count the hops to the null parent, skipping the rest of the while loop
							newlength++;
							sChunk=sChunk.parent;
						}
						length+=newlength;
						System.out.println("Coalescence time = " + (gens-length) + " generations ago.");
						new_clade_val=clade_val.substring(0,clade_val.lastIndexOf(":")+1) + length;
						//System.out.println(new_clade_val);
						pointer_map.put(sChunk,new_clade_val);											//exit method with null value assigned
						pointer_map.remove(SelectedChunkArr.get(arr_Iter));
						break; //breaks out of the main loop and writes the tree

					}
				}
			}																															


			pChunk = SelectedChunkArr.get(arr_Iter).parent;								//get the parent of the selected chunk
			
			if(arr_Iter>(selectedArrSize-1)){
				cur_tip = pointer_map.get(SelectedChunkArr.get(arr_Iter));
			}
			else{
				int tipname_advance=2;						//this is kind of lame fix to advance the tipname from A to B when there's only 1 allele since it's incremented with arr_iter which can't change.
				int arr_Iter_mod = arr_Iter%52;															//repeat from A with an incremented number when reaching end of alphabet
				if(arr_Iter_mod==0){
					subname=String.valueOf(Integer.parseInt(subname)+1);			//BUG: this probably doesn't work anymore because of skipping tips, however there's rarely more than 26 tips chosen.
				}
				
				String suffix="_0";	
				if(allele_count==1){
					suffix="_" + allele_part;																					//alleles, 0 for chrom 0, 1 for chrom 1
					if(allele_part.equals("r")){
						suffix ="_" + Integer.toString(random_allele_choice1d[arr_Iter]);
						//System.out.println("suffix=" + suffix);
					}
				}
				else{
					if(arr_Iter%2==1){
						suffix="_1";
					}
					tipname_advance=1;
				}
				cur_tip = SelectedTipNames.get(arr_Iter_mod) + "_" + tips.substring(arr_Iter_mod*tipname_advance,arr_Iter_mod*tipname_advance+1)+subname+suffix;
			}

			String length="";

			if(pointer_map.get(pChunk)==null){														//look in the hashmap for it, no pointer found
				if(cur_tip.indexOf(",")==-1){													 			//if a comma does not exist in the selectedChunk value
					int index_of_length=cur_tip.lastIndexOf(":")+1;	

					//length = cur_tip.substring(index_of_length);
					//System.out.println("DEBUG1 " + length);

					if(index_of_length==0){																		//if there's no branch length assigned yet
						length=":1";
						pointer_map.put(pChunk,cur_tip+length); 			   				//if pointer is not there, add the current tip to the hashmap
					}
					else{																											//if a branch length exists increment it by one
						length=String.valueOf((Integer.parseInt(cur_tip.substring(index_of_length)))+1);
						pointer_map.put(pChunk,cur_tip.substring(0,index_of_length)+length); 				//if pointer is not there, add the current tip to the hashmap
					}					
				}
				else{																													//comma exists

					//the methods noLength() and isCommaBetweenClades() run only if cur_tip.substring(0,1).equals("(")
					if(!cur_tip.substring(0,1).equals("(") || isCommaBetweenClades(cur_tip)){// || noLength(cur_tip) ){	//noLength() may be redundant now that commaBetweenClades is run... needs further testing
						//System.out.println(cur_tip);
						pointer_map.put(pChunk,"("+cur_tip+"):1");
						//System.out.println("("+cur_tip+"):1");
					}
					else{																												//find the final branch length and add 1
						int index_of_length=cur_tip.lastIndexOf(":")+1;
						length=String.valueOf((Integer.parseInt(cur_tip.substring(index_of_length)))+1);
						pointer_map.put(pChunk,cur_tip.substring(0,index_of_length)+length);
					}

				}
				SelectedChunkArr.add(pChunk);																//and put its pointer in the selectedChunkArr
																									 									//so it will be looked at in the next gen
			}																				
			else{ 																												//parent already found in hash map, so update it with new value
				clade_val = pointer_map.get(pChunk);
				if(cur_tip.indexOf(",")==-1){						 										//comma does not exist in selected chunk value. 
					int index_of_length=cur_tip.lastIndexOf(":")+1;

					//length = cur_tip.substring(index_of_length);
					//System.out.println("DEBUG2 " + length);

	
					if(index_of_length==0){																		//if there's no branch length assigned yet

						length=":1";
						new_clade_val = clade_val + "," + cur_tip + length;
					}	
					else{																											//if a branch length exists increment it by one
						length=String.valueOf((Integer.parseInt(cur_tip.substring(index_of_length)))+1);
						new_clade_val = clade_val + "," + cur_tip.substring(0,index_of_length)+length;
					}
				}	
				else{																												//comma exists
					//new_clade_val = clade_val + "," + "(" + cur_tip+")";		//this was the double parenthesis bug 
					new_clade_val = clade_val + ","  + cur_tip;
					//System.out.println("event");
				}	
				pointer_map.put(pChunk,new_clade_val);
			}

//UNCOMMENT TO VIEW STRING CONSTRUCTION
//			System.out.println(arr_Iter + " " + pointer_map.get(pChunk));
			pointer_map.remove(SelectedChunkArr.get(arr_Iter));						//remove current Chunk from hashmap after it's parent is added
			arr_Iter++;		

		}
		System.out.println();

		if(pointer_map.size()!=1){
			System.out.println("Taxa Did Not Coalesce!!!");
			coalesce = false;
		}
		return coalesce;
	}

	static boolean noLength(String cur_tip){
		int index_of_length=cur_tip.lastIndexOf("):")+2;						//find the last occurence of ):
		String length=cur_tip.substring(index_of_length);

		for(int i=0;i<length.length();i++){
			if(!Character.isDigit(length.charAt(i))){									//if any of the characters after ): are not digits then it needs a new length added "(...):1"
				//isNoLengthCounter++;
				return true; //if we ever see a non-digit we can end the loop early
			}
		}
		return false;
	}


	static boolean isCommaBetweenClades(String cur_tip){
		int inside_paren=0;
		for(int i=0; i<cur_tip.length(); i++){
			//iterate through cur_tip looking for a comma that's between two clades. i.e. (A:5,B:5):20,(C:5,D:5):10
			if(cur_tip.substring(i,i+1).equals("(")) { 
				inside_paren++;
			}
			else if(cur_tip.substring(i,i+1).equals(")")) { 
				inside_paren--;
			}
			else if(inside_paren==0 && cur_tip.substring(i,i+1).equals(",")){
				//isCommaCounter++;
				return true; //if we ever see a comma while inside_paren is 0 we can end the loop early
			}
		}
		return false;
	}

	static int[][] getSpeciationArray(int nInds, String foldername){

		//get this info from guide tree eventually
		GuideTree gt = new GuideTree();		

		int[][] speciationGroups = gt.getSpeciationGroups(nInds, foldername);
		//int[][] speciationGroups = {{5000, 0,6, 7,13, 14,19},{-999999}};
		//int[][] speciationGroups = {{5000, 0,9, 10,19},{7300, 0,9, 10,19, 20,39},{8000, 0,19, 20,39, 40,49, 50,59},{9000, 0,19, 20,39, 40,59, 60,66, 67,73, 74,79},{-999999}};

		//System.out.println("-------------------------------------------");
		for(int i=0; i<speciationGroups[0].length; i++){
			for(int j=0; j<speciationGroups[i].length; j++){
				//System.out.print(speciationGroups[i][j] + " ");
			}
			//System.out.println();
		}
		//System.out.println("-------------------------------------------");

		int[][] speciationArray = speciationGroups;

		return speciationArray;

	}



	static void writeToFile(String treestring, int nInds, int gens, int numNuc, double recombRateAll, int[][] speciation_time, String replicate){
	
		String speciation_time_label = "";
		int iter=0;
		while(speciation_time[iter][0]!=-999999){
			speciation_time_label+="_" + String.valueOf(gens-speciation_time[iter][0]);
			iter++;
		}

		File dir = new File("T"+speciation_time_label);
		if(!dir.exists()) {
			try{
				dir.mkdir();
				System.out.println("Directory created at: T_" +speciation_time_label );
			}
			catch(SecurityException se)
			{
				System.out.println("Directory not created");
			}
		}



		try{
			File file = new File("T" + speciation_time_label + "/tree_" + nInds + "_" + gens + "_" + numNuc + "_" + recombRateAll + "__" + replicate + ".tre");
			if(!file.exists() ){
				file.createNewFile();
			}
			BufferedWriter bw = new BufferedWriter(new FileWriter(file));
			bw.write(treestring);
			bw.close();
			System.out.println("Tree written to "+file);
		} catch(IOException e){
			e.printStackTrace();
		}
	}		

	static int[][] chooseRandomTips(int nInds, int species, int tips){	//selects a random n tips per species
		Random r = new Random();
		int[][] selectedTips = new int[species][tips];
		int[] tempArray = new int[tips];		

		for(int j=0; j<species; j++){
			for(int i=0; i<tips; i++){
				tempArray[i] = r.nextInt(nInds) + (nInds*j);
				//selectedTips[j][i]=r.nextInt(nInds) + (nInds*j);
			}
			Arrays.sort(tempArray);
			for(int i=0; i<tips; i++){
				selectedTips[j][i]=tempArray[i];
			}
		}
		return selectedTips;

	}

	static int[] splitInds(int nInds, int divisions){
	//This method splits an int into d parts evenly then add the remainder across the parts
	// e.g. 20/3 = 7 7 6
		//System.out.println("inside splitPop\n\tspecies=" + divisions + "\n");
		if(divisions==0){ int parts[] = {0, nInds}; return parts;}
		
		int n = nInds;
		int d = divisions;
		int[] parts = new int[d+1];
		//System.out.println(d);

		parts[0]=0;

		int m=n/d;
		for(int i=1; i<=d; i++){
			parts[i] = m;		
		}

		int r=n%d;
		for(int i=1;i<=r;i++){
			parts[i] += 1;
		}

		for(int i=0; i<=d; i++){
			System.out.print("\t"+parts[i]);		
		}	
		System.out.println();

		return parts;		


	}


}

class Ind{

	//int numNuc = 6;
 	//Chunk[][] Chrom = new Chunk[numNuc][2];	
	Chunk[][] Chrom;	


	Ind(int numNuc){	//root individuals
		this.Chrom = new Chunk[numNuc][2];	

		for(int j=0;j<2;j++){//two homologous chromosomes
			for(int i=0;i<numNuc;i++){
				this.Chrom[i][j]=new Chunk();
			}
		}
	}	

	Ind(Ind[] Cur_ancestors, int numNuc, double[] recombRate, Random r){
		int[][] recombArray = randomRecombination(numNuc,recombRate,r); 
		this.Chrom = new Chunk[numNuc][2];	

//uncomment to view recombination
/*
		//System.out.println(ancestors_id[0] + "," + ancestors_id[1]);
		for(int j=0;j<2;j++){
			for(int i=0;i<numNuc;i++){
				System.out.print(recombArray[i][j]);
			}
			System.out.println();
		}
*/

		Ind CurAncestor;
		for(int j=0;j<2;j++){
			CurAncestor = Cur_ancestors[j];
			for(int i=0;i<numNuc;i++){
				//Chrom[i][j]=new Chunk();
				this.Chrom[i][j]=new Chunk(CurAncestor.Chrom[i][recombArray[i][j]]);	
				//this long line creates a new chunk obj and adds it to the Chrom array. The parameter is the parent which is a Chunk in the same homologous position on that ancestor.
			}
		}
	}

	public int[][] randomRecombination(int numNuc, double[] recombRate, Random r){  //recombination is based on distance between chunks, given by recombRate array

		int[][] recombArray = new int[numNuc][2];
		//Random r = new Random();
		//System.out.println("Recomb Rate: " + recombRate);


		for(int j=0;j<2;j++){
			if(r.nextDouble()<.5){	//first allele is a coin flip
				recombArray[0][j]=0;
			}
			else{
				recombArray[0][j]=1;
			}

			if(recombRate[0]==0){
				for(int i=1;i<numNuc;i++){				
					recombArray[i][j] = recombArray[0][j];
				}
			}
			else{			
				for(int i=1;i<numNuc;i++){																			//subsequent alleles are based on recombination rate derived from distances between chunks
					if(r.nextDouble()<recombRate[i]){															//e.g. if .3 is less than .5 (random number between 0 and 1)
						recombArray[i][j] = (recombArray[i-1][j])^1;								//XOR with 1 flips the 0 to a 1 or 1 to a 0
					}
					else{
						recombArray[i][j] = recombArray[i-1][j];
					}
				}
			}
		}
		return recombArray;

	}

}
	
class Chunk{
	Chunk parent=null;

	Chunk(){
	}

	Chunk(Chunk newparent){
		this.parent=newparent;
	}
}
		
class NewThread implements Callable<ArrayList<Ind>>{

	int id;
	//int[] randomArray;
	int splitadd;
	int split_start;
	int split_end;
	int min;
	int max;
	int numNuc;
	double[] recombRate;
	private ArrayList<Ind> CurGenArr;
	ArrayList<Ind> splitArray;
	//Thread t;
	
	NewThread(int id, int split_start, int split_end, int splitadd, int min, int max, int numNuc, double[] recombRate, ArrayList<Ind> CurGenArr){
		//this.id = id;
		this.splitadd=splitadd;
		this.split_start=split_start;
		this.split_end=split_end;
		this.min = min;
		this.max = max;		
		this.numNuc = numNuc;
		this.recombRate = recombRate;
		this.CurGenArr = CurGenArr;
		this.splitArray = new ArrayList<Ind>((split_end-split_start)+1);
		//this.randomArray = randomArray;
		//t=new Thread(this, "Thread_"+id);
		//t.start();
		//System.out.println(split_end-split_start + " " + splitadd);

	}

	//ArrayList<Ind> getThreadArray(){
		//System.out.println("@@@@@@@@@@@@@@@@@@@ " + this.splitArray.size() );
	//	return splitArray;
	//}


	public ArrayList<Ind> call(){
		//System.out.println("limits " + min +" " +max);
		int[] ancestors_id = new int[2];
		//ArrayList<Ind> splitArray = new ArrayList<Ind>();
			//for the split calculate the Ind's ancestors
			//add them to the threadArray

				Random r =  new Random();
 				//long start = System.nanoTime();
				for(int i=split_start;i<=split_end;i++){
					ancestors_id[0]=randomMating(r, min,max,-0.1, splitadd); //-0.1 is so module 2 gives neither 1 or 0 in randomMating
					ancestors_id[1]=randomMating(r, min,max,ancestors_id[0], splitadd);
					//System.out.println(RecombWithTrees.myArray.length);
					//System.out.println(this.id);

					//ancestors_id[0]=randomArray[i]; //testing to see how much time is spent in randomMating()
					//ancestors_id[1]=randomArray[i+1];
					Ind[] Cur_ancestors = {CurGenArr.get(ancestors_id[0]),CurGenArr.get(ancestors_id[1])};
					//this.splitArray.add(new Ind(numNuc) );
					this.splitArray.add(new Ind(Cur_ancestors, numNuc, recombRate, r) );

				}
				//long time = System.nanoTime() - start;
		    //System.out.println("time splits " + (double)time);//(split_end-split_start));

		    //System.out.println("Thread " + id + " time total: " + (double)time);
		    //System.out.println("Thread " + id + " splits: " + (split_end-split_start+1));
			//this.splitArray = splitArray;
		return splitArray;		
	}
	
	
	static int randomMating(Random r, int min, int max, double prevInd, int splitadd){
	//if the 2nd ind selected matches the sex of the 1st individual we shift the selected individual up by one with a periodic boundary. If NUMTHREADS divides nInds into splits of odd numbers it's increased by 2 so we don't get 2 inds that cannot mate.
		int ancestor_id =	r.nextInt(max - min + 1) + min;								//min is 0 max is splitsize

		if( (ancestor_id%2) == (prevInd%2) ){														//evens are females, odds are males
			//System.out.println(ancestor_id + " " + prevInd + " " + splitadd);

			ancestor_id+=splitadd;	//add 0 or 1 rather depending on even or odd splits. 

			if(ancestor_id>max){	//if it goes beyond the max mod it by the max
				ancestor_id=(ancestor_id%max)+min;
			}

			//ancestor_id=randomMating(r, min,max,ancestor_id);
			//RecombWithTrees.myCounter++;
		}

		return ancestor_id;
	}
}
