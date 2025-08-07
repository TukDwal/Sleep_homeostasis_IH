# Sleep_homeostasis_IH
This repository makes available the code that was developed to produce the results from the paper "Dysregulation of slow wave activity in Idiopathic Hypersomnia", Tugdual Adam, Lucie Barateau, Yves Dauvilliers
Institut des Neurosciences de Montpellier, INSERM; Unité des Troubles du sommeil et de l'éveil, CHU Gui de Chauliac, Montpellier.


Raw anonymized data is not currently available as it represents nearly 1To of data. It may be partially provided upon reasonable request from any qualified investigator.
For any inquiries, please contact: tugdual.adam@gmail.com

Description:
- MAIN.m contains the main pipeline of data analysis, organized by figure number in the original paper
The other functions are supporting functions necessary for the code to run:
- mtspecgramc.m: multitaper spectral analysis
- fit_S_dual_stages.m: mathematical model of sleep pressure
- finddatagroups: supporting function to extract starting and ending indexes of consecutive values in a vector 
