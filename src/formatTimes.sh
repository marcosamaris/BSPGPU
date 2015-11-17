#!/bin/bash -xe


LOG_DIR="/home/marcos/Dropbox/Doctorate/Results/2015/BSPGPU/src/"
cd $LOG_DIR

rm *.txt    
touch matMul-Gm-SP.txt
touch matMul-Gm-Un-SP.txt
touch matMul-Sm-SP.txt        
touch matMul-Sm-Un-SP.txt        
touch SubSeqMax.txt
cd experiments
for time in $(ls -d */);
	do        
	cd $time
	for file in $(ls );
	do    
	    [ "$file" == 'matMul-Gm-SP.txt' ] && cat matMul-Gm-SP.txt | grep float | awk '{print $6}' > ../matMul-Gm-SP
	    
	    [ "$file" == 'matMul-Gm-Un-SP.txt' ] && cat matMul-Gm-Un-SP.txt | grep float | awk '{print $6}' > ../matMul-Gm-Un-SP
	    
	    [ "$file" == 'matMul-Sm-SP.txt' ] && cat matMul-Sm-SP.txt | grep float | awk '{print $6}' > ../matMul-Sm-SP
	    
	    [ "$file" == 'matMul-Sm-Un-SP.txt' ] && cat matMul-Sm-Un-SP.txt | grep float | awk '{print $6}' > ../matMul-Sm-Un-SP
	    
	    [ "$file" == 'SubSeqMax.txt' ] && cat SubSeqMax.txt | grep int | awk '{print $6}' > ../SubSeqMax
	done
cd ../        	

paste matMul-Gm-SP matMul-Gm-SP.txt > coco
mv  coco matMul-Gm-SP.txt

paste matMul-Gm-Un-SP matMul-Gm-Un-SP.txt > coco
mv  coco matMul-Gm-Un-SP.txt

paste matMul-Sm-SP matMul-Sm-SP.txt  > coco
mv coco matMul-Sm-SP.txt        

paste matMul-Sm-Un-SP matMul-Sm-Un-SP.txt > coco
mv  coco  matMul-Sm-Un-SP.txt        

paste SubSeqMax SubSeqMax.txt > coco
mv  coco SubSeqMax.txt
		
done
 	
