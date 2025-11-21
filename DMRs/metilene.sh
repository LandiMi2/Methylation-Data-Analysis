#Prepare input files for metilene
#CpG
metilene_input.pl --in1 ../../filtered/Mock/CpG/Mock2.CpG.bedGraph,../../filtered/Mock/CpG/Mock3.CpG.bedGraph \
--in2 ../../filtered/CMV/CpG/CMV2.CpG.bedGraph,../../filtered/CMV/CpG/CMV3.CpG.bedGraph --h1 mock --h2 cmv \
-o mock_cmv_cpg.input

#run metilene
metilene -a mock -b cmv mock_cmv_cpg.input | sort -V -k1,1 -k2,2n > mock_cmv_cpg.dmrs
#filter 10% methylation differences 
awk '($5 >= 10 || $5 <= -10)' mock_cmv_cpg.dmrs > mock_cmv_cpg.0.10.dmrs

#CHG
metilene_input.pl --in1 ../../filtered/Mock/CHG/Mock2.CHG.bedGraph,../../filtered/Mock/CHG/Mock3.CHG.bedGraph \
--in2 ../../filtered/CMV/CHG/CMV2.CHG.bedGraph,../../filtered/CMV/CHG/CMV3.CHG.bedGraph --h1 mock --h2 cmv \
-o mock_cmv_chg.input

metilene -a mock -b cmv mock_cmv_chg.input | sort -V -k1,1 -k2,2n > mock_cmv_chg.dmrs
awk '($5 >= 10 || $5 <= -10)' mock_cmv_chg.dmrs > mock_cmv_chg.0.10.dmrs

#CHH
metilene_input.pl --in1 ../../filtered/Mock/CHH/Mock2.CHH.bedGraph,../../filtered/Mock/CHH/Mock3.CHH.bedGraph \
--in2 ../../filtered/CMV/CHH/CMV2.CHH.bedGraph,../../filtered/CMV/CHH/CMV3.CHH.bedGraph --h1 mock --h2 cmv \
-o mock_cmv_chh.input

metilene -a mock -b cmv mock_cmv_chh.input | sort -V -k1,1 -k2,2n > mock_cmv_chh.dmrs
awk '($5 >= 10 || $5 <= -10)' mock_cmv_chh.dmrs > mock_cmv_chh.0.10.dmrs
