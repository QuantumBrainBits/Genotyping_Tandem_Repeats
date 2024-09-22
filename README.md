Genotyping Tandem Repeats

![Screenshot from 2024-09-23 00-47-11](https://github.com/user-attachments/assets/f97a206a-a574-4a01-a86e-52f51a2b3a3f)






Data Requirement:

To genotype STRs, HipSTR requires NGS data. However, as the depth of sequencing and the read length in these datasets can vary dramatically, here we briefly describe key factors to consider before generating data for HipSTR analyses.

Process of Read to Variant : stages where artifacts gets introduced.

![Screenshot from 2024-09-23 00-39-08](https://github.com/user-attachments/assets/21948e3e-09fb-4a17-b7fb-67396eb6b77f)


Because of the repetitive nature of STRs, reads that do not fully extend across the repeat only provide a lower bound on its length. While this lower bound is informative and is leveraged by HipSTR, obtaining accurate and robust STR genotypes requires reads that fully extend across the repetitive sequence (i.e. spanning reads). The number of reads that span an STR is a function of the read length, the sequencing depth, and the length of the repeat (as well as various other factors). The interplay between these factors is relatively complex, but Figure 2 in a recent review by Press et al. nicely highlights these dependencies. As one would intuitively expect, using longer read lengths and higher sequencing coverage increases the number of spanning reads. Conversely, increasing the length of the repeat reduces the number of spanning reads, making it more difficult to accurately genotype long STRs. When the number of spanning reads approaches single digits, you statistically run the risk of observing reads from only 1 out of 2 chromosome copies, making it impossible to correctly call both alleles in a heterozygous individual.


Statistical model to capture artifacts from the reads known as "Stutter error" effects more in Tandem repeat regions:

![Screenshot from 2024-09-23 00-39-23](https://github.com/user-attachments/assets/dc46a22c-8289-4746-8b4e-4bec05af9b59)

