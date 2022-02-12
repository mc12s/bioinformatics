//usage java ConcatSeqs Concat/T_10000_2500/1_12/r2000/64/u_1.5.ali 12
//use in a bash loop for each u_*.ali file within 1_12

import java.io.*;
import java.util.HashMap;
import java.util.Map;

public class ConcatSeqs{

	public static void main(String[] args){

		int seqs_per_ali = 23;

		String outFile = args[0]; // Concat/T_40000_2500/1_12/r0/64/u_1.5.ali

		String[] parts = args[0].split("/");
		int max_exon = Integer.parseInt(args[1]);

		String prefix = parts[1];
		String suffix = parts[3] + "/" + parts[4] + "/" + parts[5];

		int len = Integer.parseInt(parts[4]);

		int total_len = len * max_exon;

		HashMap<Integer, String> hmap = new HashMap<Integer, String>();

		try{
			BufferedReader[] br = new BufferedReader[12];

			for(int exon=0; exon<max_exon; exon++){
				String filename = prefix + "/" + (exon+1) + "/" + suffix;
				//System.out.println(filename);
				br[exon] = new BufferedReader(new FileReader(filename));
			}

			String full_ali = "";

			for(int ali=0; ali<100; ali++){	
				for(int seq=1; seq<=24; seq++){
					hmap.put(seq, "");	//fill hashmap to avoid checking if a value exists later
				}
				
				full_ali = full_ali + " 24 " + total_len;

				for(int exon=0; exon<max_exon; exon++){
					int count = ali*seqs_per_ali;
					String line;
					br[exon].readLine(); //skip header
					//System.out.println((exon+1) + " " +line);
					while(count <= ((ali+1)*seqs_per_ali)){
						line = br[exon].readLine();
						//System.out.println((exon+1) + " " +line);
						count++;
						String[] line_parts = line.split(" ");
						//System.out.println(line_parts[0] + " " + line_parts[line_parts.length-1]);
						hmap.put(Integer.parseInt(line_parts[0]), hmap.get(Integer.parseInt(line_parts[0])) + line_parts[line_parts.length-1]);
						//add to hashmap taxon id = key, sequence = value
						//if they appear concatenate the sequences
					}
				}
				//add these 25 lines to String
				for (int key : hmap.keySet()) {
			    full_ali += "\n" + key + " " + hmap.get(key);
				}
				full_ali += "\n";
				//full_ali += temp;
				hmap.clear();
			}

			for(int exon=0; exon<max_exon; exon++){
				br[exon].close();
			}

			//System.out.println(full_ali);
			//write full_ali to file
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outFile)));
			bw.write(full_ali);

			bw.flush();
			bw.close();
		

		}catch(IOException ioe){
			System.out.println(ioe);
		}
	}

}
