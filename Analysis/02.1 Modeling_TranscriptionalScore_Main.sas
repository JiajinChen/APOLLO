PROC IMPORT OUT= WORK.EXP 
            DATAFILE= "D:\APOLLO\TCGA_main_FPKM.txt"
            DBMS=TAB REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

PROC PHREG data=exp;
	model OS*Censor(0)= AFF1 ARHGEF10 BCR BRCA1 BRCA2 CASP8 CDK6 CHEK2 CHIC2 CNTRL EML4 ERBB2 FANCD2 IGF2BP2 ITGAV LATS2 MSH6 MSN NFE2L2 PBRM1 PLCG1 PTPN13 REL RNF213 SEPT9 SMO TNC WRN
                         / selection=stepwise slentry=0.05001
                           slstay=0.05 details; 
run;
/* CHIC2 IGF2BP2 ITGAV MSN PLCG1 */