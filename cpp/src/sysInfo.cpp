#include "../include/sysInfo.h"
#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <sys/types.h>
#include <ctime>
#include <chrono>
#include <sys/utsname.h>

void runApp(const char *appname)
{
    char *vmName = getenv("HOSTNAME");
    char *vmUserName = getenv("USER");
    std::cout << "This is the " << appname << " project running on " << vmName << " under supervision of usr: " << vmUserName << "\n";
}

void getOSDetails()
{
    char hostName[512];
    std::cout << "UUID: " << getuid();
    gethostname(hostName, sizeof(hostName) - 1);
    std::cout << hostName;
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
}

//shows the end time of the application
//Works with chrono
void showDate()
{
    auto finishTime = std::chrono::system_clock::now();
    time_t showDate = std::chrono::system_clock::to_time_t(finishTime);
    std::cout << ctime(&showDate) << "\n";
}