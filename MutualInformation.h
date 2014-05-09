// Utility Functions that can be called from a digital micrograph script inside the program.

float MutualInformation( DM_ImageToken image_token, DM_ImageToken image_token2 , long a , long b);

float MutualInformation2( DM_ImageToken image_token, DM_ImageToken image_token2 , long a , long b);

float MutualInformationCL( DM_ImageToken image_token, DM_ImageToken image_token2 , long a , long b);

void MutualInformationMap( DM_ImageToken image_token, DM_ImageToken image_token2 , long a , long b);

void XCFPCFPCPCF(DM_ImageToken image_token, DM_ImageToken image_token2 , long df);

void Magnification(DM_ImageToken image_token, DM_ImageToken image_token2);

void MagnificationTest(DM_ImageToken image_token, DM_ImageToken image_token2, float expectedDF, long numberoftrials, float minscale, float maxscale);