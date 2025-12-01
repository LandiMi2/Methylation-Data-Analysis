# Inputs are from pre-filtered coverage files (cov > 5 sites) and NOT common sites files. 

#Prepare input files for metilene
#CpG
metilene_input.pl --in1 ../../filtered/CaMV/CpG/CaMV2.CpG.bedGraph,../../filtered/CaMV/CpG/CaMV3.CpG.bedGraph \
--in2 ../../filtered/Mock/CpG/Mock2.CpG.bedGraph,../../filtered/Mock/CpG/Mock3.CpG.bedGraph --h1 camv --h2 mock \
-o camv_mock_cpg.input

#sort
metilene -c 2 -a camv -b mock camv_mock_cpg.input | sort -V -k1,1 -k2,2n > camv_mock_cpg.dmrs

#filter
metilene_output.pl -q camv_mock_cpg.dmrs -o camv_mock_cpg -a camv -b mock

##CHG
metilene_input.pl --in1 ../../filtered/CaMV/CHG/CaMV2.CHG.bedGraph,../../filtered/CaMV/CHG/CaMV3.CHG.bedGraph \
--in2 ../../filtered/Mock/CHG/Mock2.CHG.bedGraph,../../filtered/Mock/CHG/Mock3.CHG.bedGraph --h1 camv --h2 mock \
-o camv_mock_chg.input

#sort
metilene -c 2 -a camv -b mock camv_mock_chg.input | sort -V -k1,1 -k2,2n > camv_mock_chg.dmrs

#filter
metilene_output.pl -q camv_mock_chg.dmrs -o camv_mock_chg -a camv -b mock


##CHH
metilene_input.pl --in1 ../../filtered/CaMV/CHH/CaMV2.CHH.bedGraph,../../filtered/CaMV/CHH/CaMV3.CHH.bedGraph \
--in2 ../../filtered/Mock/CHH/Mock2.CHH.bedGraph,../../filtered/Mock/CHH/Mock3.CHH.bedGraph --h1 camv --h2 mock \
-o camv_mock_chh.input

#sort
metilene -c 2 -a camv -b mock camv_mock_chh.input | sort -V -k1,1 -k2,2n > camv_mock_chh.dmrs

#filter
metilene_output.pl -q camv_mock_chh.dmrs -o camv_mock_chh -a camv -b mock
