INSTALLATION
------------
FARTSSIE is a standalone Excel macro-enabled file (.xlsm). Macros need to be enabled for the spreadsheet to function.
If you are downloading FARTSSIE on a Windows 11 operating system, you will need to do the following:
-Right-click on the downloaded file
-Select "Properties"
-In the "Security" section, click "Unblock", then click "Apply".
-When you open the file in Excel, click "Enable Content".

Hidden under the visual basic code of FARTSSIE is the version history:<p>

```
VERSION HISTORY
---------------
Version 3.0: Corrected routine to calculate total CV from published CI for Bioequivalence, Parallel (PowerTOST method).
Version 2.9: Two worksheets were added to estimate the combined dropout rates and buffer the sample size accordingly: "Dropouts - 2 rates", and "Dropouts - 3 rates". Added "RefScaledBE" Worksheet - sample size table for HVDs and NTIDs.
Version 2.8: FARTSSIE moves to Github! TotalCV calculation added to Parallel BE sheet. Prompt for CV on Parallel BE updated to Total CV.
Version 2.7: PowerTOST suggested example code corrected.
There were a few "quiet" updates for version 2.8 that didn't warrant a version number increase. For instance I think I finally found and corrected the Error 13 bug on some Mac platforms.
Version 2.6: Added graphs to Parallel BE worksheet, and Mac had some macro issues with Double variable type (changed to Variant).
Version 2.5: Removed Reference-Scaled BE Limits (skewed results at higher CVs, improper method implemented). Added links to PowerTOST.
Version 2.4: Fixed the way that d was depicted on superiority and non-inferiority trials.
Version 2.3: Added Lehr's Formula to parallel superiority
Version 2.2: Added d to the first two graphs, fixed 30% threshold BE limits for HVD (limits would go within 80-125 below 30% ISV)
Version 2.1: Added snazzy distribution graphs
Version 2.0: Added TPD HVD limits, Policy on Bioequivalence Standards for HVD Products 18-Apr-16:
Version 1.9: FDA sigmaw0 changed to 0.294, EMEA widened CI limits constrained to 69.84-143.19%
Version 1.8: Fixed precision of ISV and ratio for BE sheets, updated to .xlsm format
Version 1.7: Added Reference-Scaled BE Limits to Replicate Spreadsheet for FDA
Version 1.6: Added Scaled BE Limits spreadsheet based on FDA recommendations
Version 1.5:
- added power calculator button
Version 1.4: Now compatible with Office 2007. Dr. Russel Lenth generously provided the library subroutines to calculate
- non-central distributions (NCt)
- Warning added on bioequivalence spreadsheets for ratio outside BE limits
- Limit warning at 1000 subjects to prevent needless waiting

```
  
I have also added a few random tools I have developed over the years that have been sitting on my hard drive (again, offered as a "free" second opinion).
These programs are unvalidated, only lightly tested, and come with no guarantees of cross-platform compatibility or accuracy. Some of them are being used as proof-of-concept programs in the bioequivalence courses I teach at the Leslie Dan Faculty of Pharmacy.

For instance, the bootstrapping programs I posted do not run a single ANOVA. Think "high-school" level standard error. However, they do give surprisingly good agreement with calculated power, and bootstrapping is a stochastic method after all.

“You know," said Arthur, "it's at times like this, when I'm trapped in a Vogon airlock with a man from Betelgeuse, and about to die of asphyxiation in deep space that I really wish I'd listened to what my mother told me when I was young."
"Why, what did she tell you?"
"I don't know, I didn't listen.”
― Douglas Adams, The Hitchhiker's Guide to the Galaxy
(also applicable to my relationship with statistics)
