Hidden under the visual basic code of FARTSSIE is the version history:

'VERSION HISTORY
'---------------
'Version 1.4: Now compatible with Office 2007. Dr. Russel Lenth generously provided the library subroutines to calculate
'- non-central distributions (NCt)
'- Warning added on bioequivalence spreadsheets for ratio outside BE limits
'- Limit warning at 1000 subjects to prevent needless waiting
'Version 1.5:
'- added power calculator button
'Version 1.6: Added Scaled BE Limits spreadsheet based on FDA recommendations
'Verstion 1.7: Added Reference-Scaled BE Limits to Replicate Spreadsheet for FDA
'Version 1.8: Fixed precision of ISV and ratio for BE sheets, updated to .xlsm format
'Version 1.9: FDA sigmaw0 changed to 0.294, EMEA widened CI limits constrained to 69.84-143.19%
'Version 2.0: Added TPD HVD limits, Policy on Bioequivalence Standards for HVD Products 18-Apr-16:
'Version 2.1: Added snazzy distribution graphs
'Version 2.2: Added d to the first two graphs, fixed 30% threshold BE limits for HVD (limits would go within 80-125 below 30% ISV)
'Version 2.3: Added Lehr's Formula to parallel superiority
'Version 2.4: Fixed the way that d was depicted on superiority and non-inferiority trials.
'Version 2.5: Removed Reference-Scaled BE Limits (skewed results at higher CVs, improper method implemented). Added links to PowerTOST.
'Version 2.6: Added graphs to Parallel BE worksheet, and Mac had some macro issues with Double variable type (changed to Variant).
'Version 2.7: PowerTOST suggested example code corrected.

I have also added a few random tools I have developed over the years that have been sitting on my hard drive (again, offered as a "free" second opinion).
These programs are unvalidated, only lightly tested, and come with no guarantees of cross-platform compatibility or accuracy. Some of them are being used as proof-of-concept programs in the bioequivalence courses I teach at the Leslie Dan Faculty of Pharmacy.

For instance, the bootstrapping programs I posted do not run a single ANOVA. Think "high-school" level standard error. However, they do give surprisingly good agreement with calculated power, and bootstrapping is a stochastic method after all.
