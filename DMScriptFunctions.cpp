#include "stdafx.h"
#include "DMScriptFunctions.h"


DigitalMicrograph::ScriptObject DMScript::GetFrameSetInfo (DigitalMicrograph::ScriptObject &Acquis)
{
	static DigitalMicrograph::Function __sFunction = (DM_FunctionToken) NULL;
	static const char *__sSignature = "ScriptObject GetFrameSetInfo( ScriptObject )";

	Gatan::PlugIn::DM_Variant params[2];

	params[1].v_object = (DM_ObjectToken) Acquis.get();
	GatanPlugIn::gDigitalMicrographInterface.CallFunction( __sFunction.get_ptr(), 2, params, __sSignature );

	return (DM_ScriptObjectToken_1Ref) params[0].v_object;
};

DigitalMicrograph::ScriptObject DMScript::GetFrameSetInfoPtr (DigitalMicrograph::ScriptObject &Acquis)
{
	static DigitalMicrograph::Function __sFunction = (DM_FunctionToken) NULL;
	static const char *__sSignature = "ScriptObject GetFrameSetInfoPtr( ScriptObject * )";

	Gatan::PlugIn::DM_Variant params[2];

	params[1].v_object_ref = (DM_ObjectToken*) Acquis.get_ptr();
	GatanPlugIn::gDigitalMicrographInterface.CallFunction( __sFunction.get_ptr(), 2, params, __sSignature );

	return (DM_ScriptObjectToken_1Ref) params[0].v_object;
}

void DMScript::EMGetStageXY (double &StageX, double &StageY)
{
	static DigitalMicrograph::Function __sFunction = (DM_FunctionToken) NULL;
	static const char *__sSignature = "void EMGetStageXY( double *, double* )";

	Gatan::PlugIn::DM_Variant params[2];

	params[0].v_float64_ref =  &StageX;
	params[1].v_float64_ref =  &StageY;
	GatanPlugIn::gDigitalMicrographInterface.CallFunction( __sFunction.get_ptr(), 2, params, __sSignature );
}

void DMScript::EMSetStageXY (double StageX, double StageY)
{
	static DigitalMicrograph::Function __sFunction = (DM_FunctionToken) NULL;
	static const char *__sSignature = "void EMSetStageXY( double , double )";

	Gatan::PlugIn::DM_Variant params[2];

	params[0].v_float64 =  StageX;
	params[1].v_float64 =  StageY;
	GatanPlugIn::gDigitalMicrographInterface.CallFunction( __sFunction.get_ptr(), 2, params, __sSignature );
}

void DMScript::EMSetStageX (double StageX)
{
	static DigitalMicrograph::Function __sFunction = (DM_FunctionToken) NULL;
	static const char *__sSignature = "void EMSetStageX( double )";

	Gatan::PlugIn::DM_Variant params[2];

	params[0].v_float64 =  StageX;
	GatanPlugIn::gDigitalMicrographInterface.CallFunction( __sFunction.get_ptr(), 1, params, __sSignature );
}

void DMScript::EMSetStageY (double StageY)
{
	static DigitalMicrograph::Function __sFunction = (DM_FunctionToken) NULL;
	static const char *__sSignature = "void EMSetStageY( double )";

	Gatan::PlugIn::DM_Variant params[2];

	params[0].v_float64 =  StageY;
	GatanPlugIn::gDigitalMicrographInterface.CallFunction( __sFunction.get_ptr(), 1, params, __sSignature );
}

double DMScript::EMGetStageX ()
{
	static DigitalMicrograph::Function __sFunction = (DM_FunctionToken) NULL;
	static const char *__sSignature = "double EMGetStageX( )";

	Gatan::PlugIn::DM_Variant params[1];
	GatanPlugIn::gDigitalMicrographInterface.CallFunction( __sFunction.get_ptr(), 1, params, __sSignature );
	return params[0].v_float64;
}

double DMScript::EMGetStageY ()
{
	static DigitalMicrograph::Function __sFunction = (DM_FunctionToken) NULL;
	static const char *__sSignature = "double EMGetStageY( )";

	Gatan::PlugIn::DM_Variant params[1];
	GatanPlugIn::gDigitalMicrographInterface.CallFunction( __sFunction.get_ptr(), 1, params, __sSignature );
	return 	params[0].v_float64;
}

void DMScript::EMSetHighTensionOffset( double offset )
{
    static Gatan::DM::Function __sFunction = (DM_FunctionToken) NULL;
    static const char *__sSignature = "void EMSetHighTensionOffset( double )";

    Gatan::PlugIn::DM_Variant params[1];

    params[0].v_float64 = offset;
    GatanPlugIn::gDigitalMicrographInterface.CallFunction( __sFunction.get_ptr(), 1, params, __sSignature );
}

double DMScript::EMGetHighTensionOffset( )
{
    static Gatan::DM::Function __sFunction = (DM_FunctionToken) NULL;
    static const char *__sSignature = "double EMGetHighTensionOffset(  )";

    Gatan::PlugIn::DM_Variant params[1];
    GatanPlugIn::gDigitalMicrographInterface.CallFunction( __sFunction.get_ptr(), 1, params, __sSignature );

    return params[0].v_float64;
}