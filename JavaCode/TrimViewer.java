import java.awt.*;
import java.awt.event.*;
import java.io.*;
import javax.swing.*;
import javax.swing.event.*;
import java.util.Arrays;

public class TrimViewer extends JPanel implements MouseListener{

	public TrimViewer(){
  addMouseListener(this); 

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
						br.close();
					}catch(IOException ioe){
						System.out.println("<<!!ERROR main()!!>> MESSAGE:"+ioe.getMessage());
					}				
				}
			}
		}


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
			System.out.println(Arrays.toString(maxTally));		
			}
		}
		for (int[] element : smallestInd) {
   		System.out.println(Arrays.toString(element));
		}		

		for (double[] element : windowTally) {
   		System.out.println(Arrays.toString(element));
		}		


		JCheckBox textCheckbox = new JCheckBox("Show Text");
    textCheckbox.setSelected(false);
		textCheckbox.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e) {         
			        Object source = e.getItemSelectable();
				if(e.getStateChange() == ItemEvent.DESELECTED){
					showtext=false;
					repaint();
				}
				else{				
					showtext=true;				
					repaint();				
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
					repaint();
				}
				else{				
					prop_param=0;
					mingoodsites_param=0;
					untrimmed=true;				
					repaint();				
				}
			}
		});


		JPanel jplCheckBox = new JPanel();
		jplCheckBox.setLayout(new BorderLayout() );
			setBorder(BorderFactory.createEmptyBorder(490,-800,0,0));

		jplCheckBox.add(textCheckbox);
		add(jplCheckBox);

		JPanel jplCheckBox2 = new JPanel();
		jplCheckBox2.setLayout(new BorderLayout() );
			setBorder(BorderFactory.createEmptyBorder(490,-550,0,0));

		jplCheckBox2.add(untrimmedCheckbox);
		add(jplCheckBox2);

		JSlider zoomSlider = new JSlider(JSlider.HORIZONTAL, 0, 10, 5);
		
		zoomSlider.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent ce) {
				baseSize=1+2*( (JSlider)ce.getSource()).getValue();
				repaint();
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
		g.setColor(Color.white);	//background for alignment			
		g.fillRect(1,490,1050,1600); 	//background for alignment
		g.setColor(Color.black);
		g.drawRect(1,490,1050,1600); 	//border for alignment

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
		

		if(untrimmed){
				seqString=untrimmed_seqString;
			}
			else{
				seqString=trimmed_seqString;
			}
			


		int	locus_start=locus_id;
		int	locus_end=locus_id2;

		if(!(locus_id<=locus_id2)){ 
			locus_start=locus_id2;
			locus_end=locus_id;		
		}

//		int seqLength=seqString[0][locus_start][prop_param][mingoodsites_param].length();
//		for(int k=locus_start+1; k<=locus_end; k++){
//			seqLength=seqLength+seqString[0][k][prop_param][mingoodsites_param].length();
//		}

		if(!loci1_selected){

			int startpos=0;
			int startrow=0;
			textwrap=0;
			for(int k=locus_start; k<=locus_end; k++){

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
								case "A":if(numA+numR>.5*ntaxa){g.setColor(Color.lightGray);}else{g.setColor(Color.blue);}break;
								case "T":if(numT+numY>.5*ntaxa){g.setColor(Color.lightGray);}else{g.setColor(Color.red);}break;		
								case "C":if(numC+numY>.5*ntaxa){g.setColor(Color.lightGray);}else{g.setColor(Color.green);}break;
								case "G":if(numG+numR>.5*ntaxa){g.setColor(Color.lightGray);}else{g.setColor(Color.yellow);}break;
								case "R":g.setColor(Color.black);break;
								case "Y":g.setColor(Color.black);break;
								default:g.setColor(Color.black);break;
							}		
						}
						else{
								g.setColor(Color.lightGray);
						}			

						//if(j*baseSize>1000){textwrap=j*baseSize/1000;}
						if((10+(baseSize*(j+startpos)))>1000){textwrap=((startpos+j)*baseSize)/1000;}

						if(untrimmed){mingoodsites_param=0;prop_param=0;}
						//System.out.println((10+(baseSize*(j-(textwrap))))+(510+(textwrap/100)+(baseSize*i)));	
						if(Character.toString(seqString[i][k][prop_param][mingoodsites_param].charAt(j)).equals("-")){ //|| Character.toString(seqString[i][locus_id].charAt(j)).equals("N") || Character.toString(seqString[i][locus_id].charAt(j)).equals("n")){
							g.setColor(Color.black);
							g.fillRect(10+(baseSize*(j+startpos))-(1000*textwrap), 525+(2*ntaxa*textwrap*baseSize)+((baseSize+1)*i),baseSize,1);		
						}
						else if(Character.toString(seqString[i][k][prop_param][mingoodsites_param].charAt(j)).equals("N")){ //|| Character.toString(seqString[i][locus_id].charAt(j)).equals("N") || Character.toString(seqString[i][locus_id].charAt(j)).equals("n")){
							g.setColor(Color.black);
							g.fillRect(10+(baseSize*(j+startpos))-(1000*textwrap), 520+(2*ntaxa*textwrap*baseSize)+((baseSize+1)*i),baseSize,baseSize);		
						}
						else{
							g.fillRect(10+(baseSize*(j+startpos))-(1000*textwrap), 520+(2*ntaxa*textwrap*baseSize)+((baseSize+1)*i),baseSize,baseSize);		
							if(showtext && baseSize > 8){					
								g.setColor(Color.black);					
								g.setFont(g.getFont().deriveFont((float)baseSize));
								g.drawString(Character.toString(seqString[i][k][prop_param][mingoodsites_param].charAt(j)),10+(baseSize*(j+startpos))-(1000*textwrap), 530+(2*ntaxa*textwrap*baseSize)+((baseSize+1)*i));	
							}
						}
					}								
				}
				startpos=startpos+seqLength; //Continue to next locus at this pixel coefficient
			}
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
	


	double largestMask=0.0;
	double largestWindow;
	double smallestWindow;
	int prop_iter=5;
	int mingoodsites_iter=9;
	int nloci=801; //CHANGE
	int ntaxa=5;		//CHANGE
	int[][] smallestInd=new int[5][9];			
	double[] maxTally=new double[ntaxa+1];			
	double[][] windowTally=new double[5][9];
	double[][][] locimasked=new double[5][9][nloci];
	double[][] smallestMaxTally=new double[5][9];
	int pointX=0;
	int pointY=0;
	int prop_param=0;
	int mingoodsites_param=0;
	int locus_id;
	int locus_id2;
	int new_locus_id;
	boolean loci1_selected=false;
	int locisqrt;
	int horiz_side;
	int vert_side;
	int squareSz;
	int lociX;
	int lociY;
	int lociX2;
	int lociY2;
	int new_lociX;
	int new_lociY;
	int baseSize=10;
	int textwrap=0;
	int currbasepos=0;
	int currbaserow=0;
	boolean showtext=false;
	boolean untrimmed=false;
	String[][][][] seqString = new String[ntaxa][nloci][5][9];
	String[][][][] trimmed_seqString = new String[ntaxa][nloci][5][9];
	String[][][][] untrimmed_seqString = new String[ntaxa][nloci][1][1];

	public void mouseEntered(MouseEvent mouse){ }   
	public void mouseExited(MouseEvent mouse){ }
	public void mousePressed(MouseEvent mouse){ }
	public void mouseReleased(MouseEvent mouse){ }

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
			repaint();			
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
					repaint();
				}		
				else{
					lociY2=new_lociY;
					lociX2=new_lociX;	
					locus_id2=new_locus_id;
					loci1_selected=false;
					repaint();	
				}
			}
		}			
	}


  public static void main(String[] a) { //Window setup 
    JFrame window = new JFrame();
    window.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    window.setBounds(30, 30, 1280, 1024);
    window.add(new TrimViewer());
    window.setVisible(true);
  }	
}
