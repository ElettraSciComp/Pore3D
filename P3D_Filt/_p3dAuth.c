#include <stdio.h>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <signal.h>
#include <netdb.h>
#include <unistd.h>
#include <string.h>
#include <sys/un.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/ioctl.h>
#include <linux/if.h>

#include <time.h>

#include "p3dAuth.h"

#define ADDRESS "dynamite.esce.elettra.trieste.it"
//#define ADDRESS "localhost"
#define PORT 5000

#define BUFFLEN 2048

static char recv_gh[BUFFLEN];

char* authenticate_real(char* funname) {
    int sock, len;
    char msg[BUFFLEN];
    char msgtmp[BUFFLEN];

    struct hostent *he;
    struct sockaddr_in serv_name;
    struct in_addr h_addr;

    int i, ct, rc;
    time_t mytime;
    struct tm *mytm;

    // clear buffer:
    memset(recv_gh, 0, BUFFLEN * sizeof (char));

    // Socket creation:
    sock = socket(AF_INET, SOCK_STREAM, 0);

    if (sock < 0) {
        close(sock);
        exit(1);
    }

    bzero(&serv_name, sizeof (serv_name));

    serv_name.sin_family = AF_INET;

    // With a direct IP address:
    //cli_name.sin_addr.s_addr = inet_addr(ADDRESS);

    // With a DNS name:
    if ((he = gethostbyname(ADDRESS)) == NULL) {
        sprintf(msgtmp, "error resolving hostname %s", ADDRESS);
        return msgtmp;
    }

    // The host might be connected to multiple networks and have different
    // addresses on each one. We scan h_addr_list (the vector is terminated
    // by a null pointer).
    ct = 0;
    while (he->h_addr_list[ct] != NULL) {
        memcpy(&serv_name.sin_addr, he->h_addr_list[ct], he->h_length);
        serv_name.sin_port = htons(PORT);

        rc = connect(sock, (struct sockaddr *) &serv_name, sizeof (serv_name));
        if (rc != -1)
            break;
        else
            ct++;
    }

    // After attempting all the IP addresses, if there is a connection error exit:
    if (rc == -1) {
        h_addr.s_addr = *((unsigned long *) he->h_addr_list[ct - 1]);
        sprintf(msgtmp, "error connecting to %s", inet_ntoa(h_addr));
        return msgtmp;
    }

    //
    // We are connected to the server. Let's send information:
    //

    // Clear buffer:
    memset(msg, 0, BUFFLEN * sizeof (char));
    len = sizeof (serv_name);

    // Getting the client MAC id for sending it via socket:
    struct ifreq s;
    int fd = socket(PF_INET, SOCK_DGRAM, IPPROTO_IP);

    strcpy(s.ifr_name, "eth0");
    if (0 == ioctl(fd, SIOCGIFHWADDR, &s)) {
        for (i = 0; i < 6; ++i) {
            //printf("%02x", (unsigned char) s.ifr_addr.sa_data[i]);
            sprintf(msgtmp, "%02x", (unsigned char) s.ifr_addr.sa_data[i]);
            strcat(msg, msgtmp);
        }
    }

    // Get the current USER for sending it via socket: (comma separated):
    sprintf(msgtmp, ",%s", getenv("USER"));
    strcat(msg, msgtmp);

    // Add the name of the FUNCTION (get from input parameter):
    sprintf(msgtmp, ",%s", funname);
    strcat(msg, msgtmp);

    // Add a TIMESTAMP:
    mytime = time(NULL);
    mytm = localtime(&mytime);

    strftime(msgtmp, sizeof msgtmp, ",%x-%X", mytm);
    strcat(msg, msgtmp);

    // Send the string via socket:
    send(sock, &msg, strlen(msg), 0);
    rc = recvfrom(sock, &recv_gh, BUFFLEN, 0, (struct sockaddr*) &serv_name, &len);
    //message received //printf("message received (lib): (%d chars) %s\n",rc,recv_gh);

    // Close socket and return:
    close(sock);
    return recv_gh;
}

char authenticate_local(char* funname)
{       
    FILE* fvol;
    char  *mac;
    char  ret_s;
    
    // Set the number of threads as half the CPU power
    // in order to avoid saturation:
    omp_set_num_threads(omp_get_num_procs()/2);
    
    // Read the MAC address on eth0:    
    if ((fvol = fopen("/sys/class/net/eth0/address", "rb")) == NULL) {
        ret_s = '0';
    }
    // Read raw data from file:
    mac = (char*) malloc(17*sizeof(char));
    fread(mac, sizeof (char), 17, fvol);
    
    // Perform the comparison:
    if (strcmp(mac,"00:23:ae:fe:ed:d4") == 0)
        ret_s = '1';
    else
        ret_s = '0';    
    
    /* Close file handler: */
    fclose(fvol);
    free(mac);    
    
    return ret_s;
}

char authenticate(char* funname)
{              
    return '1';
}