#include "hystreApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

// App revision
#include "HystreRevision.h"

InputParameters
hystreApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy DirichletBC, that is, set DirichletBC default for preset = true
  params.set<bool>("use_legacy_dirichlet_bc") = false;

  return params;
}

hystreApp::hystreApp(InputParameters parameters) : MooseApp(parameters)
{
  hystreApp::registerAll(_factory, _action_factory, _syntax);
}

hystreApp::~hystreApp() {}

void
hystreApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"hystreApp"});
  Registry::registerActionsTo(af, {"hystreApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
hystreApp::registerApps()
{
  registerApp(hystreApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
hystreApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  hystreApp::registerAll(f, af, s);
}
extern "C" void
hystreApp__registerApps()
{
  hystreApp::registerApps();
}
