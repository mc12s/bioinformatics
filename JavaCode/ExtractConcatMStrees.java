//This program will extract the ms recombination trees and output them to a new file for further processing
//usage java ExtractConcatMStrees T_40000_2500/10000reps_0_12.msout 1 6

import java.io.*;

public class ExtractConcatMStrees{

	public static void main(String[] args){

		String inFile = args[0];
		int exonA = Integer.parseInt(args[1]);
		int exonB = Integer.parseInt(args[2]);
		String newfilename = inFile + "extracted" + exonA + "_" + exonB + ".txt";


		try{
		
			File file_in = new File(inFile);
			if(!file_in.exists()){
				System.out.println("File not found. Exiting...");
				System.exit(0);
			}

			File file = new File(newfilename);
			if (!file.exists()) {
		     file.createNewFile();
			}
			BufferedWriter bw = new BufferedWriter(new FileWriter(file));
			BufferedReader br = new BufferedReader(new FileReader(inFile));	

			br.readLine();
			br.readLine();
			br.readLine();	//skip header
			String line = br.readLine();
			String prevline = "";

			int len = 0;
			boolean noskip = false;
			boolean skipA = false;

			int repcount = 0;
			//int skippedreps = 0;
			//startrep = startrep - 100;

			while(line != null){
				if(line.equals("//")){	//start of replicate
					/*if(skippedreps < startrep){
						//System.out.println("here");
						line = br.readLine();				
						skippedreps++;
						continue;
					}*/
					len=0;
					noskip = true;
					skipA = false;
					line = br.readLine();	//get first line

					//System.out.println(line);
					//bw.write(line);
					//bw.newLine();
				}
				//System.out.println("test: " + skippedreps);
				if(noskip){
					int start = line.indexOf('[');
					int end = line.indexOf(']');
					if(start!=-1){
						//System.out.println(line);
						len+=Integer.parseInt(line.substring(start+1,end));
						//System.out.println(len);
						if(len>=exonA && !skipA){
							//repcount++;	
							//System.out.println(line);
							bw.write(line.split("]")[1]);
							bw.newLine();
							skipA = true;						
						}

						if(len>=exonB){
							repcount++;	
							//System.out.println("in " + line);
							bw.write(line.split("]")[1]);
							bw.newLine();
							noskip = false;
						}
					}
					if(line.equals("")){	//reached end of replicate without finding brackets (only 1 tree)
						repcount++;
						bw.write(prevline);	//exonA
						bw.newLine();
						bw.write(prevline); //exonB same tree owing to zero recombination
						bw.newLine();
					}
				}
				
				prevline = line;			//store line - 1 
				line = br.readLine();	//advance

			}

			if(repcount<10000){
				bw.write(prevline);	//exonA
				bw.newLine();
				bw.write(prevline); //exonB same tree owing to zero recombination
				bw.newLine();
			}

			bw.flush();
			bw.close();
			
			System.out.println("File written to: " + newfilename);
		}catch(IOException ioe){
			System.out.println(ioe);
		}

	}

}

