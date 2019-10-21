#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <sys/types.h>

using namespace std;

void runApp(const char *appname)
{
    char *vmName = getenv("HOSTNAME");
    char *vmUserName = getenv("USER");
    cout << "This is the " << appname << " project running on " << vmName << " under supervision of usr: " << vmUserName << "\n";
}

void getOSDetails()
{
    char hostName[512];
    cout << "UUID: " << getuid();
    gethostname(hostName, sizeof(hostName) - 1);
    cout << hostName;
}

int main()

{
    const char *appname = "Wobbling fit 2 params";
    runApp(appname);
    return 0;
}