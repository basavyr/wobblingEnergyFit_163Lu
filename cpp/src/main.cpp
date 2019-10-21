#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <sys/types.h>
#include <ctime>
#include <chrono>
#include <sys/utsname.h>

//optional OMP
#include<omp.h>

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

void sysInfo()
{
    struct utsname buffer;

    errno = 0;
    if (uname(&buffer) != 0)
    {
        perror("uname");
        exit(EXIT_FAILURE);
    }

    printf("system name = %s\n", buffer.sysname);
    printf("node name   = %s\n", buffer.nodename);
    printf("release     = %s\n", buffer.release);
    printf("version     = %s\n", buffer.version);
    printf("machine     = %s\n", buffer.machine);

#ifdef _GNU_SOURCE
    printf("domain name = %s\n", buffer.domainname);
#endif

    // return EXIT_SUCCESS;
}

void runParallel()
{
#pragma omp parallel for
for(int i=0;i<4;++i)
cout<<i<<" "<<i+1<<"\n";
#pragma omp parallel

  {
    printf("https://helloacm.com\n");
  }
}



int main()

{
    const char *appname = "Wobbling fit 2 params";
    runApp(appname);
    auto finishTime = chrono::system_clock::now();
    time_t showDate = chrono::system_clock::to_time_t(finishTime);
    sysInfo();
    cout << ctime(&showDate) << "\n";
    runParallel(); 
   return 0;
}
