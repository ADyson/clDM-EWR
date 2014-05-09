namespace DMScript
{
	DigitalMicrograph::ScriptObject GetFrameSetInfo (DigitalMicrograph::ScriptObject &Acquis);
	DigitalMicrograph::ScriptObject GetFrameSetInfoPtr (DigitalMicrograph::ScriptObject &Acquis);
	void EMGetStageXY (double &StageX, double &StageY);
	void EMSetStageXY (double StageX, double StageY);
	void EMSetStageX (double StageX);
	void EMSetStageY (double StageY);
	double EMGetStageX ();
	double EMGetStageY ();
	void EMSetHighTensionOffset( double offset );
	double EMGetHighTensionOffset( );
}