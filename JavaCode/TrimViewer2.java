import java.awt.*;
import java.awt.event.*;
import java.io.*;
import javax.swing.*;
import javax.swing.JMenuBar;
import javax.swing.event.*;
import java.util.Arrays;


////////TRY MOVING ALIGNMENT STRINGS INTO A CHAR ARRAY TO GET TEXT TO RENDER FASTER??


public class TrimViewer2 {

	public static int locus_start=0;
	public static int locus_end=0;
	public static int baseSize=10;
	public static int textwrap=0;
	public static int currbasepos=0;
	public static int currbaserow=0;
	public static boolean showtext=false;
	public static boolean untrimmed=false;
	public static int nloci=801; //CHANGE
	public static int ntaxa=5;		//CHANGE
	public static String[][][][] seqString = new String[ntaxa][nloci][5][9];


	private double largestMask=0.0;
	private double largestWindow;
	private double smallestWindow;
	private int prop_iter=5;
	private int mingoodsites_iter=9;

	private int[][] smallestInd=new int[5][9];			
	private double[] maxTally=new double[ntaxa+1];			
	private double[][] windowTally=new double[5][9];
	private double[][][] locimasked=new double[5][9][nloci];
	private double[][] smallestMaxTally=new double[5][9];
	private int pointX=0;
	private int pointY=0;
	private int prop_param=0;
	private int mingoodsites_param=0;
	private int locus_id;
	private int locus_id2;
	private int new_locus_id;
	private boolean loci1_selected=false;
	private int locisqrt;
	private int horiz_side;
	private int vert_side;
	private int squareSz;
	private int lociX;
	private int lociY;
	private int lociX2;
	private int lociY2;
	private int new_lociX;
	private int new_lociY;

	private String[][][][] trimmed_seqString = new String[ntaxa][nloci][5][9];
	private String[][][][] untrimmed_seqString = new String[ntaxa][nloci][1][1];




	//PROCEED WITH GLOBALS FOR NOW
	public TrimViewer2(){};

	private class HeatMaps extends JComponent implements MouseListener{


		private HeatMaps(){

  	addMouseListener(this); 

		/////////////////////////////////////////////////////////////////////////////////////////////
		//HEATMAPS
		/////////////////////////////////////////////////////////////////////////////////////////////

		for (int i=0;i<prop_iter;i++){
			for(int j=0;j<mingoodsites_iter;j++){
				try{				//Load masked loci data into locimasked array to create loci masking heatmap
					BufferedReader brMissing = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File("../Results/test/P0072_"+i+"_"+j+"_lociMasked.txt") ) ));			
					for(int k=0;k<nloci;k++){
						String tempS=brMissing.readLine();
						locimasked[i][j][k]=Double.parseDouble(tempS);
						if (locimasked[i][j][k]>largestMask){largestMask=locimasked[i][j][k];}					
					}
					brMissing.close();
				}
				catch(IOException e){
					System.out.println("ERROR: No File found!");					
				}			
			}
		}
		

		smallestWindow=9999;
		largestWindow=0;
		for (int i=0;i<prop_iter;i++){
			for(int j=0;j<mingoodsites_iter;j++){				
				smallestMaxTally[i][j]=99999;
				smallestInd[i][j]=0;
				try{ 			//Load MaxTally and WindowTally results into arrays to be used for visualization
					BufferedReader brMaxTally = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File("../Results/P0072_"+i+"_"+j+"_maxTally.txt") ) ));			
					BufferedReader brWindowTally = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File("../Results/P0072_"+i+"_"+j+"_windowTally.txt") ) ));			

					Arrays.fill(maxTally,0);
					for (int k=0;k<nloci;k++){
						String tempS=brMaxTally.readLine();
						String[] parts=tempS.split("\t");
						for (int m=0;m<=ntaxa;m++){						
							maxTally[m]=maxTally[m]+Double.parseDouble(parts[m]);
						}
					}

					for (int k=0;k<7+j;k++){  //get the mingoodsites parameter
						brWindowTally.readLine();
					}					
					String tempS=brWindowTally.readLine();
					String[] parts=tempS.split("\t");
					windowTally[i][j]=0;
					for (int m=0;m<parts.length;m++){
						windowTally[i][j]=windowTally[i][j]+Double.parseDouble(parts[m]); //sum the mingoodsites value of all taxa
					}
					windowTally[i][j]=Math.log10(windowTally[i][j]);
					if(windowTally[i][j]<smallestWindow){smallestWindow=windowTally[i][j];}
					if(windowTally[i][j]>largestWindow){largestWindow=windowTally[i][j];}

					for (int m=0;m<=ntaxa;m++){	//find index of lowest MaxTally value(find optimal minimum proportion of same base per site e.g. .5)
						maxTally[m]=Math.log10(maxTally[m]);
						if(maxTally[m]<smallestMaxTally[i][j]){
							smallestMaxTally[i][j]=maxTally[m];
							smallestInd[i][j]=m+3;  //add 1 for index starts at 0, starts at .2
							//System.out.println(smallestInd[i][j]);
						}
					}
					//System.out.println(smallestInd.length);
					//System.out.println(Arrays.toString(smallestInd[2]));
					brMaxTally.close();
					brWindowTally.close();
				}
				catch(IOException e){
					System.out.println("ERROR: No File found!");					
				}	
			//System.out.println(Arrays.toString(maxTally));		
			}
		}
		for (int[] element : smallestInd) {
   		//System.out.println(Arrays.toString(element));
		}		

		for (double[] element : windowTally) {
   		//System.out.println(Arrays.toString(element));
		}		


		JCheckBox textCheckbox = new JCheckBox("Show Text");
    textCheckbox.setSelected(false);
		textCheckbox.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e) {         
			        Object source = e.getItemSelectable();
				if(e.getStateChange() == ItemEvent.DESELECTED){
					showtext=false;
					updatePanel();
				}
				else{				
					showtext=true;
					updatePanel();		
				}
			}
		});

		JCheckBox untrimmedCheckbox = new JCheckBox("Untrimmed Alignments");
    untrimmedCheckbox.setSelected(false);
		untrimmedCheckbox.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e) {         
			        Object source = e.getItemSelectable();
				if(e.getStateChange() == ItemEvent.DESELECTED){
					untrimmed=false;
					updatePanel();
					
				//align.updatePanel();
				}
				else{				
					prop_param=0;
					mingoodsites_param=0;
					untrimmed=true;				
					updatePanel();				
				}
			}
		});


		JPanel jplCheckBox = new JPanel();
		jplCheckBox.setLayout(new FlowLayout() );
		//setBorder(BorderFactory.createEmptyBorder(490, -880,0,0));

		jplCheckBox.add(textCheckbox);
		add(jplCheckBox);

		JPanel jplCheckBox2 = new JPanel();
		jplCheckBox2.setLayout(new BorderLayout() );
			setBorder(BorderFactory.createEmptyBorder(490, 550,0,0));

		jplCheckBox2.add(untrimmedCheckbox);
		add(jplCheckBox2);

		JSlider zoomSlider = new JSlider(JSlider.HORIZONTAL, 0, 10, 5);
		
		zoomSlider.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent ce) {
				baseSize=1+2*( (JSlider)ce.getSource()).getValue();
				updatePanel();
			}
		});		
		
		zoomSlider.setMinorTickSpacing(2);
		zoomSlider.setMajorTickSpacing(10);

		zoomSlider.setLabelTable(zoomSlider.createStandardLabels(10));

		add(zoomSlider, BorderLayout.CENTER);
  }

	
	public void paintComponent(Graphics g) {
	
		super.paintComponent(g);
		//Generates loci masking heatmap - From light blue to black


		squareSz=15;
		locisqrt=(int)Math.ceil(Math.sqrt(nloci)); //determine best dimensions for heatmap
		int locicount=0;

	  double colorrange1=255/largestMask;

		g.setColor(Color.gray);
		g.fillRect(squareSz+20, squareSz, locisqrt*squareSz, locisqrt*squareSz);
		outerLoop: //paint the locimasked_map
		for(int i=0;i<locisqrt;i++){
			g.setColor(Color.black);		
			g.drawString(String.valueOf(i*locisqrt+1),2,25+15*i);
			for(int j=0;j<locisqrt;j++){
				if(locicount>nloci-1){break outerLoop;}			
				int aColor=(int)(locimasked[prop_param][mingoodsites_param][i*locisqrt+j]*colorrange1);
				//System.out.println(aColor+" "+(i*locisqrt+j));
				Color c = new Color(0, 0, 255-aColor, 255);
				g.setColor(c);
				g.fillRect (squareSz+20+j*squareSz, squareSz+i*squareSz, squareSz, squareSz);  
				locicount++;			
			}
		}

		g.setColor(Color.black);	//border for loci window
		g.drawRect(squareSz+20, squareSz, locisqrt*squareSz, locisqrt*squareSz);
		

		//setup size of parameter heatmap to scale with the size of the loci heatmap !CHECK FOR LARGE LOCI SIZE
		int boxside=squareSz*locisqrt;
		vert_side=boxside/prop_iter;
		horiz_side=boxside/mingoodsites_iter;
		double colorrange2=200/(largestWindow-smallestWindow);

		//paint the parameters heatmap
		for(int i=0;i<prop_iter;i++){
			g.setColor(Color.black);
			g.drawString("--------",squareSz+locisqrt*15+horiz_side+40,19+vert_side*i);
			g.drawString("0."+String.valueOf(i+1),squareSz+locisqrt*15+horiz_side+30,65+vert_side*i);
			for(int j=0;j<mingoodsites_iter;j++){
				int aColor=(int)((windowTally[i][j]-smallestWindow)*colorrange2);
				Color c = new Color(0+aColor, 255-aColor, 0, 255);
				g.setColor(c);
				g.fillRect (squareSz+locisqrt*15+100+j*horiz_side, 15+i*vert_side, horiz_side, vert_side); 
			}
		}

		g.setColor(Color.black);
		for(int j=0;j<mingoodsites_iter;j++){
			g.drawString(String.valueOf(j+8),squareSz+locisqrt*15+horiz_side*(j+2)+horiz_side/2,30+prop_iter*vert_side);				
		}

		if(!loci1_selected){
			g.setColor(Color.red);	
			g.drawRect(squareSz+20+lociX2*squareSz, squareSz+lociY2*squareSz, squareSz, squareSz);  				
		}

		g.setColor(Color.white);
		g.drawRect(squareSz+20+lociX*squareSz, squareSz+lociY*squareSz, squareSz, squareSz);  

		g.setColor(Color.black);
		g.drawRect(squareSz+locisqrt*15+100, 15+vert_side*(smallestInd[0][0]-1), horiz_side*mingoodsites_iter, vert_side); 
		g.drawRect(squareSz+locisqrt*15+100-1, 15+vert_side*(smallestInd[0][0]-1)-1, horiz_side*mingoodsites_iter+2, vert_side+2); 
		g.setColor(Color.white);	
		g.drawRect(squareSz+locisqrt*15+100+1+(horiz_side*mingoodsites_param), 15+(vert_side*prop_param)+1, horiz_side-1, vert_side-1); 

		g.setColor(Color.black); //border for parameters window
		g.drawRect(squareSz+locisqrt*15+100, squareSz, horiz_side*mingoodsites_iter, vert_side*prop_iter);
		
		//Setup nLoci graph, missingalowed
		g.setColor(Color.black); //border for parameters window
		g.drawRect(squareSz+locisqrt*35+10, squareSz, horiz_side*mingoodsites_iter, vert_side*prop_iter);
		g.setColor(Color.white); //border for parameters window
		g.fillRect(squareSz+locisqrt*35+11, squareSz+1, horiz_side*mingoodsites_iter-2, vert_side*prop_iter-2);


		locus_start=locus_id;
		locus_end=locus_id2;

		if(!(locus_id<=locus_id2)){ 
			locus_start=locus_id2;
			locus_end=locus_id;		
		}

		g.setColor(Color.black);		
		g.setFont(g.getFont().deriveFont((float)12));
		g.drawString("Locus id:"+Integer.toString(locus_id+1),700,510);			
		g.drawString("-",785,510);			
		
		if(!loci1_selected){
			g.setColor(Color.red);		
			g.drawString(Integer.toString(locus_id2+1),800,510);			
		}
	}
	

	public void mouseEntered(MouseEvent mouse){ }   
	public void mouseExited(MouseEvent mouse){ }
	public void mousePressed(MouseEvent mouse){ }
	public void mouseReleased(MouseEvent mouse){ }


	private void updatePanel(){
			revalidate();
			repaint();		
	}

	public void mouseClicked(MouseEvent mouse){
	 		
		pointX = mouse.getX();
		pointY = mouse.getY();
		
		//System.out.println(pointY+" "+pointX);

		if(pointX>550 && pointX<980 && pointY>16 && pointY<456){ //clicked parameters window
			if(untrimmed){			
				prop_param=0;
				mingoodsites_param=0;
			}
			else{			
				prop_param=(pointY-16)/vert_side;
				mingoodsites_param=(pointX-550)/horiz_side;
			}			
			//System.out.println(prop_param+" "+mingoodsites_param);			
			updatePanel();		//////////////////////////////////////////?FIX 
			//align.updatePanel();		
		}			

		if(pointX>35 && pointX<470 && pointY>16 && pointY<436){		//Clicked loci window
			new_lociY=(pointY-16)/squareSz;
			new_lociX=(pointX-35)/squareSz;
			new_locus_id=(new_lociY*locisqrt)+new_lociX;

			//System.out.println(locus_id);	
			if(new_locus_id<nloci){	
				if(!loci1_selected){
					lociY=new_lociY;
					lociX=new_lociX;	
					locus_id=new_locus_id;
					loci1_selected=true;
					updatePanel();	
					//align.updatePanel();
				}		
				else{
					lociY2=new_lociY;
					lociX2=new_lociX;	
					locus_id2=new_locus_id;
					loci1_selected=false;
					updatePanel();
					//align.updatePanel();
				}
			}
		}
	}
	}





	public class Alignment extends JComponent{
		//TrimViewer2 trimviewer;
		int[][][] trimmed_seqLength = new int[nloci][prop_iter][mingoodsites_iter];
		int totSeqLength;
		char seqSelected[][];

		public char[][] getSeqSelected(){
			totSeqLength=0;

			for(int k=locus_start; k<=locus_end; k++){
				totSeqLength = totSeqLength+trimmed_seqLength[k][prop_param][mingoodsites_param];
			}

			char[][] seqSelected=new char[ntaxa][totSeqLength];

			if(!loci1_selected){
				for(int k=locus_start; k<=locus_end; k++){

					int curSeqLength = trimmed_seqLength[k][prop_param][mingoodsites_param];
					int oldLength=0;
					System.out.println(curSeqLength);
					for(int i=0; i<ntaxa; i++){
						for(int j=0; j<curSeqLength; j++){
							//System.out.println("Zero");
							//System.out.println(seqString[i][k][prop_param][mingoodsites_param].charAt(j));
							seqSelected[i][j+oldLength]=seqString[i][k][prop_param][mingoodsites_param].charAt(j);
						}
					}
					oldLength=oldLength+curSeqLength;
					System.out.println(oldLength);
				}
				//System.out.println("xxxx:"+seqSelected[0][0]);
			}
			return seqSelected;
		}


		protected void updateAlignPanel(){
			//char testChar[][]= testGetChar();
			//System.out.println(testChar[0][0]);
			seqSelected = getSeqSelected();
			System.out.println("Check4");
			System.out.println(seqSelected[0][0]);

			System.out.println("Check3");
			revalidate();
			repaint();
		}

		Color alignColor[][][][][] = new Color[ntaxa][nloci][5][9][3000];//FIXME maybe just do this for the selected prop_iter...

		public Alignment(){

			//JLabel label1 = new JLabel("Image and Text");
		  //add(label1);
		  System.out.println("check");
			
			addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent e) {
				System.out.println("Alignment pane was clicked at : " + e.getPoint());
				updateAlignPanel();

				}
			});


			for(int k=0; k<nloci; k++){
				for(int i=0; i<prop_iter; i++){
					for(int j=0; j<mingoodsites_iter; j++){			
						String alignmentFile = "../TrimmedAlignments/P0072_"+i+"_"+j+"_L"+(k+1)+".fasta";
						int taxa_id=0;
						
						try{
							BufferedReader br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(alignmentFile) ) ));
							String tempS=br.readLine();
							while(tempS!=null){
								tempS=br.readLine(); //skip header
								trimmed_seqString[taxa_id][k][i][j]=tempS;		
								taxa_id++;
								tempS=br.readLine();
							}
							trimmed_seqLength[k][i][j]=trimmed_seqString[0][k][i][j].length();
							//System.out.println("Check5");
							br.close();
						}catch(IOException ioe){
							System.out.println("<<!!ERROR main()!!>> MESSAGE:"+ioe.getMessage());
						}				
					}
				}
			}
	
			//int[][][] untrimmed_seqLength;
			for(int k=0; k<nloci; k++){
				String alignmentFile = "../Alignments/P0072_L"+(k+1)+".fasta";
				int taxa_id=0;
				try{
					BufferedReader br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(alignmentFile) ) ));
					String tempS=br.readLine();
					while(tempS!=null){
						tempS=br.readLine(); //skip header
						untrimmed_seqString[taxa_id][k][0][0]=tempS;		
						taxa_id++;
						tempS=br.readLine();
					}
					br.close();
				}catch(IOException ioe){
					System.out.println("<<!!ERROR main()!!>> MESSAGE:"+ioe.getMessage());
				}				
			}

			if(untrimmed){
				seqString=untrimmed_seqString;
			}
			else{
				seqString=trimmed_seqString;
			}
		

			for(int k=0; k<nloci; k++){
								//System.out.println(k);
				int seqLength=seqString[0][k][prop_param][mingoodsites_param].length();

				for(int j=0; j<seqLength; j++){
					boolean sameBase=true;
					//boolean skipN=false;
					int numA=0;
					int numT=0;
					int numC=0;
					int numG=0;
					int numR=0; //G or A
					int numY=0; //T or C
					//int numW=0; //ADD W, M, V, S, H, D
					int numN=0;
							
					for(int i=0; i<ntaxa; i++){
						switch(Character.toString(seqString[i][k][prop_param][mingoodsites_param].charAt(j))){
							case "A":numA+=1;break;
							case "T":numT+=1;break;		
							case "C":numC+=1;break;
							case "G":numG+=1;break;						
							case "Y":numY+=1;break;
							case "R":numR+=1;break;
							default:numN+=1;break;

						}
					}						

					if(numA+numN+numR>=ntaxa|| numT+numN+numY>=ntaxa || numC+numN+numY>=ntaxa || numG+numN+numR>=ntaxa){
						sameBase=true;
					} 
					else{
						sameBase=false;
					}

					for(int i=0; i<ntaxa; i++){
						if(!sameBase){			
							switch(Character.toString(seqString[i][k][prop_param][mingoodsites_param].charAt(j))){
								case "A":if(numA+numR>.5*ntaxa){alignColor[i][k][prop_param][mingoodsites_param][j]=Color.lightGray;}else{alignColor[i][k][prop_param][mingoodsites_param][j]=Color.blue;}break;
								case "T":if(numT+numY>.5*ntaxa){alignColor[i][k][prop_param][mingoodsites_param][j]=Color.lightGray;}else{alignColor[i][k][prop_param][mingoodsites_param][j]=Color.red;}break;
								case "C":if(numC+numY>.5*ntaxa){alignColor[i][k][prop_param][mingoodsites_param][j]=Color.lightGray;}else{alignColor[i][k][prop_param][mingoodsites_param][j]=Color.green;}break;
								case "G":if(numG+numR>.5*ntaxa){alignColor[i][k][prop_param][mingoodsites_param][j]=Color.lightGray;}else{alignColor[i][k][prop_param][mingoodsites_param][j]=Color.yellow;}break;
								case "R":alignColor[i][k][prop_param][mingoodsites_param][j]=Color.black;break;
								case "Y":alignColor[i][k][prop_param][mingoodsites_param][j]=Color.black;break;
								default:alignColor[i][k][prop_param][mingoodsites_param][j]=Color.black;break;
							}		
						}
						else{
								alignColor[i][k][prop_param][mingoodsites_param][j]=Color.lightGray;
						//		untrimmed_alignColor[i][k][prop_param][mingoodsites_param][j]=Color.lightGray;
						}			
					}
				}
			}
			seqSelected = getSeqSelected();
		}
		//char[][] seqSelected=getSeqSelected();
		
		public void paintComponent(Graphics g) {
			super.paintComponent(g);
			baseSize=10;
			System.out.println("CHECK2");
		  int visibleBound[]=getVisibleBound();
			int viewport_x=visibleBound[0];
			int viewport_width=visibleBound[1];
			int viewport_end=viewport_x+viewport_width;

			g.setColor(Color.white);	//background for alignment			
			g.fillRect(1,10,totSeqLength*baseSize,300); 	//background for alignment
			g.setColor(Color.black);
			g.drawRect(1,9,totSeqLength*baseSize,302); 	//border for alignment
			g.drawString(Integer.toString(locus_start),1,10);
			g.drawString(Integer.toString(locus_end),200,10);




		/////////////////////////////////////////////////////////////////////////
		//////////////////////PAINT THE SEQUENCE
		/////////////////////////////////////////////////////////////////////////
			if(!loci1_selected){
				//char seqSelected[][];// = getSeqSelected();
				int startpos=0;
				int startrow=0;
				textwrap=0;
				//for(int k=locus_start; k<=locus_end; k++){
				System.out.println("Total:" +totSeqLength);
				System.out.println("Seq:" +seqSelected[1][1]);
				
				for(int i=0; i<ntaxa; i++){
					for(int j=0; j<totSeqLength; j++){
						//int seqLength=seqString[0][k][prop_param][mingoodsites_param].length();
						if(j*baseSize>=viewport_x-25 && j*baseSize<viewport_end){
						//for(int j=0; j<seqLength; j++){

							g.setColor(Color.blue);
							//FIXME g.setColor(alignColor[i][k][prop_param][mingoodsites_param][j]);
							g.fillRect(10+(baseSize*(j+startpos)), 25+((baseSize+1)*i),baseSize,1);					

							//if(j*baseSize>1000){textwrap=j*baseSize/1000;}
							//if((10+(baseSize*(j+startpos)))>1000){textwrap=((startpos+j)*baseSize)/1000;}

							//if((10+(baseSize*(j+startpos)))>1000){textwrap=0;}

							if(untrimmed){mingoodsites_param=0;prop_param=0;}
							//System.out.println((10+(baseSize*(j-(textwrap))))+(510+(textwrap/100)+(baseSize*i)));	
							//if(Character.toString(seqString[i][k][prop_param][mingoodsites_param].charAt(j)).equals("-")){ //|| Character.toString(seqString[i][locus_id].charAt(j)).equals("N") || Character.toString(seqString[i][locus_id].charAt(j)).equals("n")){
							if(seqSelected[i][j] == '-'){ 
								g.setColor(Color.black);
								g.fillRect(10+(baseSize*(j+startpos)), 25+((baseSize+1)*i),baseSize,1);		
							}
							//else if(Character.toString(seqString[i][k][prop_param][mingoodsites_param].charAt(j)).equals("N")){ //|| Character.toString(seqString[i][locus_id].charAt(j)).equals("N") || Character.toString(seqString[i][locus_id].charAt(j)).equals("n")){
							else if(seqSelected[i][j] == 'N'){ //|| Character.toString(seqString[i][locus_id].charAt(j)).equals("N") || Character.toString(seqString[i][locus_id].charAt(j)).equals("n")){
								g.setColor(Color.black);
								g.fillRect(10+(baseSize*(j+startpos)), 20+((baseSize+1)*i),baseSize,baseSize);		
							}
							else{
								g.fillRect(10+(baseSize*(j+startpos)), 20+((baseSize+1)*i),baseSize,baseSize);		
								if(!showtext && baseSize > 8){					
									g.setColor(Color.black);					
									g.setFont(g.getFont().deriveFont((float)baseSize));
									g.drawChars(seqSelected[i],j,1,10+(baseSize*(j+startpos)), 30+((baseSize+1)*i));	//Character.toString(seqString[i][k][prop_param][mingoodsites_param].charAt(j)),10+(baseSize*(j+startpos))-(1000*textwrap), 30+(2*ntaxa*textwrap*baseSize)+((baseSize+1)*i));	
								}
							}
						}								
					}
					//startpos=startpos+seqLength; //Continue to next locus at this pixel coefficient
				}
				System.out.println("Final pos:" + startpos);
				//area.width=seqLength*baseSize;				
				
				setPreferredSize(new Dimension(50+totSeqLength*baseSize,400));
				revalidate();				
			}
		}
	}  



	public static void main(String[] args) {
        SwingUtilities.invokeLater(new Runnable() {

            @Override
            public void run() {
                new TrimViewer2().createAndShowUI();
            }
        });
    }

	private void createAndShowUI() {

        JFrame frame = new JFrame("Trim Viewer");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		 		JMenuBar menuBar = new JMenuBar();
				frame.setJMenuBar(menuBar);
        JMenu fileMenu = new JMenu("File");
        JMenu editMenu = new JMenu("Edit");
        menuBar.add(fileMenu);
        menuBar.add(editMenu);

        initComponents(frame.getContentPane());

        frame.pack();
        frame.setVisible(true);
    }

	public int[] getVisibleBound() {			//find the current viewport location and width
				Rectangle visibleRect = alignPanel.getVisibleRect();
				return new int [] {visibleRect.x, visibleRect.width};
	}
					
				JComponent topPanel = new HeatMaps();
				JComponent alignPanel = new Alignment();

	public void initComponents(Container contentPane) {

				JPanel mainPanel = new JPanel(new GridLayout(), false);

	
				JScrollPane topScrollPanel = new JScrollPane(topPanel);
				JScrollPane alignScrollPanel = new JScrollPane(alignPanel);

        topPanel.setPreferredSize(new Dimension(1536, 520));
        alignPanel.setPreferredSize(new Dimension(1536, 520));

        JSplitPane jsp = new JSplitPane(JSplitPane.VERTICAL_SPLIT);

        jsp.add(topScrollPanel, JSplitPane.TOP);
        jsp.add(alignScrollPanel, JSplitPane.BOTTOM);

        mainPanel.add(jsp);
        contentPane.add(mainPanel);

    }
}
