//This code reads alignments generated by seqgen and counts the unique sequences
//Usage: java UniqueSeqs u_1.5.ali

import java.util.*;
import java.io.*;

public class UniqueSeqs{

	public static void main(String[] args){

		//BufferedReader for ali files
		String inFile = args[0];
		//System.out.println("UniqueSeqs\tDistance");

		try{
			BufferedReader br = new BufferedReader(new FileReader(inFile));	
			int count = 0;
			int total_distance = 0;

			br.readLine();	//skip header
			String line = br.readLine();

			Map<String,String> seqMap = new HashMap<String,String>();

			while(line != null){
				//try to add to hashmap
				if(line.substring(0,3).equals(" 24")){
					//System.out.println("Unique seqs= " + seqMap.size());
					count = count + seqMap.size();

					total_distance += CalculateDistances(seqMap);

					seqMap.clear();
					line = br.readLine();
					continue;
				}

				String parts[] = line.split(" ");
				
				line = parts[parts.length-1];
				if(seqMap.get(line)==null){
					seqMap.put(line,line);
				}
				//System.out.println(line);	
				line = br.readLine();
			}
		
			//Final Replicate			
			count = count + seqMap.size();
			total_distance += CalculateDistances(seqMap);

			System.out.print(inFile+"\t");
			//System.out.println("Unique seqs= " + seqMap.size());
			//System.out.println("Total unique seqs = " + count);
			System.out.print((((double)count)/100.0));
			//System.out.println(((double)count)/100.0);
			System.out.print("\t\t");
			System.out.println((((double)total_distance)/100.0));

			seqMap.clear();
			
			

		}catch(IOException ioe){
			System.out.println(ioe);
		}

	}

	public static int CalculateDistances(Map<String,String> seqMap){
		//for each entry in the hashmap calculate the hamming distance to all others.
		int distance = 0;
		String[] values = seqMap.values().toArray(new String[seqMap.size()]);
		for (int i = 0; i < values.length; i++) {	//each value
			for (int j = i+1; j<values.length; j++) {	// is compared to the others
				for (int k = 0; k < values[i].length(); k++){	// for the length of the string
			    		if (values[i].charAt(k) != values[j].charAt(k)) {	// +1 distance if a char is different
		      				distance++;
					}
		    		}
		  	}
		}

		return distance;

	}

}
