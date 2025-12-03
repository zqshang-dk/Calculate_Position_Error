//char strEpheNam0[] = {"E:\\temp99\\P0201880.18N"};
char strEpheNam1[] = { "E:\\STUDY\\Sophomore1\\GNSS\\work3\\SJ111708S.23P" };
int EphemerisBlockNum;
EPHEMERISBLOCK *pGpsEphemeris=NULL;

ReadNavigationFile(strEpheNam1, EphemerisBlockNum, pGpsEphemeris);


//基本广播星历块
struct EPHEMERISBLOCK
//每小时一个卫星对应一个基本星历块
{
	//PRN号 
	int PRN;
	//星历参考时间
	int year, month, day, hour, minute;
	double second;
	double a0,a1,a2;//时间改正数
	//六个轨道参数
	double IODE,Crs,Deltan,M0;// ORBIT - 1
	double Cuc,e,Cus,SqrtA;// ORBIT - 2
	double Toe,Cic,OMEGA,Cis;// ORBIT - 3
	double i0,Crc,omega,OMEGAdot;// ORBIT - 4
	double IDOT,GpsWeekNumber,L2C,L2P;// ORBIT - 5
	double SatAccuracy,SatHealth,TGD,IODC;// ORBIT - 6
};

struct GPSTIME
{
	int weekno;
	double weekSecond;
};

int Calendar2GpsTime (int nYear, int nMounth, int nDay, int nHour, int nMinute, double dSecond, double &WeekSecond);
int ReadBrodcastEphemeris30(FILE *pfEph, int &EphemerisBlockNum, EPHEMERISBLOCK *pGpsEphemeris, int Flag);


//读广播星历文件，数据存储与上面定义的指针中
//参数：pfEph 表示广播星历文件指针
//      pGpsEphemeris 指向广播星历的数据指针
//      EphemerisBlockNum 返回读取到的GPS星历块个数
//      Flag=1,返回GPS星历块数，Flag=0，读取数据
int ReadBrodcastEphemeris30(FILE *pfEph, int &EphemerisBlockNum, EPHEMERISBLOCK *pGpsEphemeris, int Flag)
{
	int i, j, P1, P2, P3, P4;
	char cTemp, cSatType;
	int HeadLineNum = 1;
	int WeekNo;
	double WeekSecond;
	//定义读取的参数
	int mPrn;//卫星号PRNo 
	int year, month, day, hour, minute;//卫星钟参考时刻
	double second;
	double a0, a1, a2;//卫星钟飘参数
	double IODE, Crs, DeltN, M0;//数据星历发布时间，在轨道径向方向上周期改正正弦的振幅
	double Cuc, e, Cus, sqrtA;//轨道延迹方向上周期改正余弦振幅 、扁心率、轨道延迹方向上周期改正正弦振幅 、长半轴平方根 
	double Toe, Cic, OMEGA, Cis;//星历参考时刻、轨道倾角周期改正余弦项振幅、参考时刻升交点赤径主项、轨道倾角周期改正正弦项振幅
	double i0, Crc, omega, OMEGADOT;//参考时间轨道倾角、在轨道径向方向上周期改正余余弦的振幅、近地点角距、升交点赤径在赤道平面中的长期变化
	double IDOT, L2C, GPSWeek, L2P;////轨道倾角变化率、？？、gps周
	double AccuracyofSat, HealthofSat, TGD, IODC;//卫星精度、卫星健康、电离层群迟改正数
	double SendTime, Temp1, Temp2, Temp3;


	//读入头文件
	char strLine[255];
	while (!feof(pfEph))
	{
		fgets(strLine, 254, pfEph);
		HeadLineNum++;
		if (NULL != strstr(strLine, "END OF HEADER"))
			break;
	}

	if (Flag)
	{
		//计算星历块数
		int AllNum = 0;
		while (!feof(pfEph))
		{
			fgets(strLine, 254, pfEph);
			sscanf(strLine, "%c", &cSatType);
			switch (cSatType)
			{
			case 'G':
				for (i = 0; i < 7; i++)
					fgets(strLine, 254, pfEph);
				AllNum++;
				break;
			case 'R':
				for (i = 0; i < 3; i++)
					fgets(strLine, 254, pfEph);
				break;
			case 'C':
				for (i = 0; i < 7; i++)
					fgets(strLine, 254, pfEph);
				break;
			case 'E':
				for (i = 0; i < 7; i++)
					fgets(strLine, 254, pfEph);
				break;
			default:
				break;
			}
		}
		//临时读入星历块
		EphemerisBlockNum = AllNum;
		return 0;
	}

	//pGpsEphemeris = new EPHEMERISBLOCK[EphemerisBlockNum];
	//GPSTIME  *pGpsTime = new GPSTIME[EphemerisBlockNum];
	//fseek(pfEph, 0, SEEK_SET);
	//while (!feof(pfEph))
	//{
	//	fgets(strLine, 254, pfEph);
	//	if (NULL != strstr(strLine, "END OF HEADER"))
	//		break;
	//}

	i = 0;
	while (!feof(pfEph))
	{
		fgets(strLine, 254, pfEph);
		sscanf(strLine, "%c", &cSatType);
		switch (cSatType)
		{
		case 'G':
		{
			//读取卫星PRN号，星历参考时间
			sscanf(strLine, "%c %d %d %d %d %d %d %lf %lf %lf %lf", &cTemp, &mPrn, &year, &month, &day, &hour, &minute, &second, &a0, &a1, &a2);

			pGpsEphemeris[i].PRN = mPrn;
			pGpsEphemeris[i].year = year;
			pGpsEphemeris[i].month = month;
			pGpsEphemeris[i].day = day;
			pGpsEphemeris[i].hour = hour;
			pGpsEphemeris[i].minute = minute;
			pGpsEphemeris[i].second = second;

			//WeekNo = Calendar2GpsTime(year, month, day, hour, minute, msecond, WeekSecond);
			//pGpsTime[i].weekno = WeekNo;
			//pGpsTime[i].weekSecond = WeekSecond;
			//读取卫星钟差改正参数
			//fscanf(pfEph, "%lf %lf %lf", &a0, &a1, &a2);
			
			pGpsEphemeris[i].a0 = a0;
			pGpsEphemeris[i].a1 = a1;
			pGpsEphemeris[i].a2 = a2;
			//读 IODE,Crs,DeltN,M0
			fscanf(pfEph, "%lf %lf %lf %lf", &IODE, &Crs, &DeltN, &M0);
			pGpsEphemeris[i].IODE = IODE;
			pGpsEphemeris[i].Crs = Crs;
			pGpsEphemeris[i].Deltan = DeltN;
			pGpsEphemeris[i].M0 = M0;

			//读 Cuc,e,Cus,sqrtA
			fscanf(pfEph, "%lf %lf %lf %lf", &Cuc, &e, &Cus, &sqrtA);
			pGpsEphemeris[i].Cuc = Cuc;
			pGpsEphemeris[i].e = e;
			pGpsEphemeris[i].Cus = Cus;
			pGpsEphemeris[i].SqrtA = sqrtA;
			//Toe,Cic,OMEGA,Cis;
			fscanf(pfEph, "%lf %lf %lf %lf", &Toe, &Cic, &OMEGA, &Cis);
			pGpsEphemeris[i].Toe = Toe;
			pGpsEphemeris[i].Cic = Cic;
			pGpsEphemeris[i].OMEGA = OMEGA;
			pGpsEphemeris[i].Cis = Cis;
			//i0,Crc,w,OMEGADOT
			fscanf(pfEph, "%lf %lf %lf %lf", &i0, &Crc, &omega, &OMEGADOT);
			pGpsEphemeris[i].i0 = i0;
			pGpsEphemeris[i].Crc = Crc;
			pGpsEphemeris[i].omega = omega;
			pGpsEphemeris[i].OMEGAdot = OMEGADOT;

			//IDOT,L2Cod,GPSWeek,L2PCod
			fscanf(pfEph, "%lf %lf %lf %lf", &IDOT, &L2C, &GPSWeek, &L2P);
			pGpsEphemeris[i].IDOT = IDOT;
			pGpsEphemeris[i].L2C = L2C;
			pGpsEphemeris[i].L2P = L2P;
			pGpsEphemeris[i].GpsWeekNumber = GPSWeek;
			//AccuracyofSat,HealthofSat,TGD,IODC
			fscanf(pfEph, "%lf %lf %lf %lf", &AccuracyofSat, &HealthofSat, &TGD, &IODC);
			pGpsEphemeris[i].SatAccuracy = AccuracyofSat;
			pGpsEphemeris[i].SatHealth = HealthofSat;
			pGpsEphemeris[i].TGD = TGD;
			pGpsEphemeris[i].IODC = IODC;
			
			fscanf(pfEph, "%lf %lf %lf %lf", &SendTime, &Temp1, &Temp2, &Temp3);
			i++;
		}
			break;
		case 'R':
			for (j = 0; j < 3; j++)
				fgets(strLine, 254, pfEph);
			break;
		case 'C':
			for (j = 0; j < 7; j++)
				fgets(strLine, 254, pfEph);
			break;
		case 'E':
			for (j = 0; j < 7; j++)
				fgets(strLine, 254, pfEph);
			break;
		default:
			break;
		}
	}

	return 0;
}



//从年 月 日 转换为GPS week and week second
int Calendar2GpsTime(int nYear, int nMounth, int nDay, int nHour, int nMinute, double dSecond, double &WeekSecond)
{
	int4 DayofMonth = 0;
	int4 DayofYear = 0;
	int4 weekno = 0;
	int4 dayofweek;
	int4 m;
	if (nYear < 1980 || nMounth < 1 || nMounth > 12 || nDay < 1 || nDay > 31)  return weekno;
	//计算从1980年到当前的前一年的天数
	for( m = 1980 ; m < nYear ; m++ )
	{
		if ( (m%4 == 0 && m%100 != 0) || (m%400 == 0) ) 
		{
			DayofYear += 366;
		}
		else
			DayofYear += 365;
	}
	//计算当前一年内从元月到当前前一月的天数
	for( m = 1;m < nMounth; m++)
	{
		if(m==1 || m==3 || m==5 || m==7 || m==8 || m==10 || m==12)
			DayofMonth += 31;
		else if (m==4 || m==6 || m==9 || m==11) 
			DayofMonth += 30;
		else if (m ==2)
		{
			if ( (nYear%4 == 0 && nYear%100 != 0) || (nYear%400 == 0) )
				DayofMonth += 29;
			else 
				DayofMonth += 28;
				
		}
	}
	DayofMonth = DayofMonth + nDay - 6;//加上当月天数/减去1980年元月的6日
		
	weekno = (DayofYear + DayofMonth) / 7;//计算GPS周
	dayofweek = (DayofYear + DayofMonth) % 7;
	//m_nDayofWeek = dayofweek;
	//计算GPS 周秒时间
	WeekSecond = dayofweek * 86400 + nHour * 3600 + nMinute * 60 + dSecond;
	
	return weekno;
}