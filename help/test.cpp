//求解社交网络的最小的正影响支配集(Positive Influence Dominating Set)的两个文献中已有的贪心近似算法、几个自行设计的新的贪心算法
//***Algorithm 1: 已有的算法1(Greedy1_PIDS): 一个简单的近似贪心算法---近似比为H(δ),这里H()是和谐函数,δ是图的最大点度; 时间复杂度为O(n^3)
//其主要策略:每次选择能提升更多个未满足点的满足度的一个非支配点变为支配点.
//来自论文: On positive influence dominating sets in social networks,Theoretical Computer Science 412 (2011) 265–269
//***Algorithm 2: 已有的算法2(Greedy2_PIDS): 另一个简单的贪心算法---没有证明近似比, 但据说效果要比上面的贪心法要好; 时间复杂度仍为O(n^3)
//其主要策略:每次选择能提升其所有邻点不满足程度和最大的一个非支配点变为支配点.
//来自论文: A New Algorithm for Positive Influence Dominating Set in Social Networks.Proc.ASONAM 2012,253-257,IEEE

//***Algorithm 3: 自行新设计的算法1(NewGreedy1_PIDS): 上述贪心算法的策略的结合(实验表明基本上能保持解的质量跟已有贪心法差不多); 时间复杂度为O(n^3)
//***Algorithm 4: 自行新设计的算法2(NewGreedy2_PIDS): 上述贪心算法的策略的另一种形式的结合(实验表明基本上能保持解的质量跟已有贪心法差不多); 时间复杂度为O(n^3)
//***Algorithm 5: 自行新设计的算法3(NewGreedy3_PIDS): 快速算法(实验证明解的质量比已有贪心法的要差一点,但运行时间快, 适用于大的社交网络); 时间复杂度为O(n^2)

//2016-04-16
//注:为了简单起见,下文中提到支配集DS(Dominating Set)都是指正影响支配集PIDS(Positive Influence Dominating Set)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <conio.h>
#include <time.h>
#include <windows.h>
#include <math.h>

#define   GMaxVertexNum 37000             //图中点数的上界---根据实际情况可调整

//***注:为了提高算法的效率, 我们对读取或者产生的图进了如下预处理: 按照节点度排序,然后按度序列由大到小重新给所有节点编号,使得新的1号节点度最大,新的最后节点度最小.
//***这样各个算法中所得的正影响支配集是节点重新编号后的结果,没有转换为原图中节的序号(也没有必要转换过去,因为实验中我们仅需要所求得的正影响支配集的大小和运行时间).


/*鲍志强
这样的话，我们只能求得正影响集大小，其他的求不了。因为每次选择一个点进入D的话，都会进行每个顶点度数的排序。
*/
bool      GAdjMatrix[GMaxVertexNum+1][GMaxVertexNum+1];//图的邻接矩阵(0-1矩阵),表示图G.为全局数组能用于所有子程序中
short int DegreeList[GMaxVertexNum+1];                 //存储图中各点的度,适合最大点度不超过2^15=32768.为全局数组能用于所有子程序中
bool      DSList[GMaxVertexNum+1];                     //存储图中被选中的正影响支配集, 第i个元素=1表示点i为支配点,=0表示点i非支配点.为全局数组能用于所有子程序中

//概念1: 每点i的满足度Satisfied[i]——i邻域中支配点个数与(d(i)+1)/2之差. 值>=0表示点i为已满足点, <0表示点i为未满足点且其绝对值大小表示未满足程度.适合最大点度不超过2^15=32768.为全局数组能用于所有子程序中
short int Satisfied[GMaxVertexNum+1];
//概念2: 每个非支配点i(即DSList[i]=0)的覆盖面Coverage[i]——i邻域中未满足点的个数.值非负,表示该点如果变为支配点则可改善这些邻点的满足度.适合最大点度不超过2^15=32768.为全局数组能用于所有子程序中
short int Coverage[GMaxVertexNum+1];      //注: 一旦非支配点i变成了支配点后(DSList[i]=1),值Coverage[i]就无意义.
//概念3: 每个非支配点i(即DSList[i]=0)的被需要程度Needed[i]——i邻域中未满足点的满足度之和(为负,其绝对值意味着点i的这些邻点需要点i变为支配点的强烈程度).为全局数组能用于所有子程序中
int Needed[GMaxVertexNum+1];              //注: 一旦非支配点i变成了支配点后(DSList[i]=1),值Needed[i]就无意义. 另,Needed[..]中有些元素值可能很大,最大可能会达到2*Delta*Delta, 其中Delta为图的最大度.

int    gvertexnum;                        //图中点的实际个数
long   gedgenum;                          //图中边的实际条数
short int maxdegree;                      //保存图的最大点度,适合最大点度不超过2^15=32768
short int mindegree;                      //图的最小度
short int evendegreenum;                  //图的度为偶数的点的个数
short int onedegreenum;                   //度为1的点的个数

//***贪心策略(包括文献中已有的两个贪心算法和自行设计的几种新的贪心算法中所采用的几种基本策略)***

//策略1: 优先选择覆盖面大(见概念2: Coverage[])的非支配点变为支配点. 含义: 能提升最多个未满足点的满足度
//策略2: 优先选择被需要程度最大(见概念3: Needed[])的非支配点变为支配点. 含义: 能缓解需要程度最高(呼声最高)的那些未满足邻点的满足度.
//策略3: 优先选择未满足程度最大的点u(见概念1:Satisfied[]),再采用适当策略(比如策略1或策略2)来选取其一个非支配点型邻点w变为支配点, 来改善u的未满足度.
//策略4: 优先选择未满足程度最大的点u(见概念1:Satisfied[]),再采用适当策略(比如策略1或策略2)来选取其若干个非支配点型邻点变为支配点, 让u点一次性成为已满足点.


//Part0: 输入输出
void ReadGraph(char *txtgraphfilename);    //从文件中读取图的数据，得到图的邻接矩阵及点度列表.
                                           //***注意: 为了能在算法中提高效率,还在该过程中按照点度由大到小重新给点编号,使得节点1的度最大,节点gvertexnum的度最小
void CreateGraph(char *txtgraphfilename);  //随机产生一个图(给定点数和平均度),并保存该图到指定的文件中
                                           //***注意: 为了能在算法中提高效率, 还在该过程中按照点度由大到小重新给点编号,使得节点1的度最大,节点gvertexnum的度最小


//Part1: ***文献中已有的两种算法***
//贪心近似算法1(文献中方法)
int Greedy1_PIDS(void);  //贪心策略: 策略1
//贪心算法2(文献中方法)
int Greedy2_PIDS(void);  //贪心策略: 策略2


//Part2: ***自行设计的新算法***
//新设计的贪心算法1
int NewGreedy1_PIDS(void); //文献中两贪心策略的结合: 先采用策略1, 在策略1遇到多个点可选的情况下在这些点中采用策略2来作选择
//新设计的贪心算法2
int NewGreedy2_PIDS(void); //文献中两贪心策略的结合: 先采用策略2, 在策略2遇到多个点可选的情况下在这些点中采用策略1来作选择
//新设计的贪心算法3
int NewGreedy3_PIDS(void); //文献中两贪心策略的结合: 先采用策略2, 在策略2遇到多个点可选的情况下在这些点中采用策略1(换为覆盖面小的优先)来作选择


//新设计的局部贪心方法1
int LocalGreedy1_PIDS(void);     //策略: 每次选择当前最不满足点u, 然后在u的非支配点型邻点中反复使用策略1直到u满足为止
//新设计的局部贪心方法2
int LocalGreedy2_PIDS(void);     //策略: 每次选择当前最不满足点u, 然后在u的非支配点型邻点中反复使用策略2直到u满足为止
//新设计的局部贪心方法3(快速算法--该算法时间复杂度为O(n^2),其他算法时间复杂都是O(n^3))
int LocalGreedy3_PIDS(void);     //策略: 每次选择当前最不满足点u, 然后在u的非支配点型邻点中按点度由大到小依次添加新的支配点直到u满足为止.(读图或建图时预处理中图点序已按点度降序排了,故自然次序即可)


//新设计的新的局部贪心方法1
int NewLocalGreedy1_PIDS(void);  //策略: 每次选择当前最不满足点u, 然后在u的非支配点型邻点中使用策略1(使用1次)
//新设计的新的局部贪心方法2
int NewLocalGreedy2_PIDS(void);  //策略: 每次选择当前最不满足点u, 然后在u的非支配点型邻点中使用策略2(使用1次)
                                 //最不满足点的非支配点型邻点中使用策略2(含义:A最需要别人+B也最被人需要)
//新设计的新的局部贪心方法3(快速算法--该算法时间复杂度为O(n^2),其它算法时间复杂都是O(n^3))
int NewLocalGreedy3_PIDS(void);  //策略: 每次选择当前最不满足点u, 然后在u的非支配点型邻点中添加最度的点为新支配点.(读图或建图时预处理中图点序已按点度降序排了)



//Part4: ***可用于各算法中的子程序***
int RefinePIDS(int dssize);    //可对任何算法得到的正影响支配集DSList[1...gvertexnum]进行极小化处理(其大小为dssize),去掉不必要的支配点最终得到更小的正影响支配集DSList[1...gvertexnum],并返回新的大小
bool CheckPIDS(int dssize);    //***用于调试阶段: 对任何算法得到的正影响支配集DSList[1...gvertexnum](其大小为dssize)进行正确性检查, 正确返回1;否则返回0.

/*
//命令行执行部分(begin)*****************************************************************************************************
void main(int argc,char *argv[])
{//命令行执行该程序, 有三个参数: 第一个参数是可执行的文件名, 第二个参数是已有的图文件名, 第三个参数是保存计算结果的文件

    if(argc!=3)
	{
       printf("The number of the parameters is wrong!\n");
	   return;
	}

	char *graphfilename;
	char *resultfilename;

	graphfilename=argv[1];    //获取: 图文件名
	resultfilename=argv[2];   //获取: 结果文件名
//命令行执行部分(end)*********************************************************************************************************

*/




//***文献中已有的贪心近似算法1***
int Greedy1_PIDS(void)
{//策略:优先选择能提升更多个未满足点的满足度(即覆盖面最大)的一个非支配点变为支配点.
 //实现方法: 优先选取数组Coverage[...]中取最大值的那个非支配节点(DSList[]值0的点是非支配点)变为支配点.
	int i,j,k,dssize;

	for(i=1;i<=gvertexnum;i++)
	{//初始化
		DSList[i]=0;                         //***记录支配点(值为0的点为非支配点,值为1的点为支配点)
		Satisfied[i]=-(DegreeList[i]+1)/2;   //***记录每个点的未满足程度(值>=0表示是已满足点; 值<0表示未满足点,且其绝对值表示其未满足的程度)
		//满意度初始化为这点度数的-1/2
		Coverage[i]=DegreeList[i];           //***记录每个非支配点(即DSList[.]值为0的点)的邻域内包含有的未满足的点个数--即记录每个非支配点的覆盖面
	}
    dssize=0;//支配点个数初始化


	//预处理(略): 所有点度1的点邻点必定选作支配点(某些度为1可能也要被选作支配点,比如一个点与多个度1的点相邻时)
	//贪心近似算法
	int maxcoverage, bestvertex;
	int uncoverednum;
	uncoverednum=gvertexnum;   //初始化,表示开始时未满足点的个数

	while(uncoverednum)        //如果还有未满足的节点就需要继续添加新的支配点
	{
		maxcoverage=0;        //初始化.该变量表示新添加一个支配点最多可支配多少个未满足的点
		for(i=1;i<=gvertexnum;i++)
			if(DSList[i]==0 && Coverage[i]>maxcoverage)  //所有非支配点都可能要选取(包括某些度为1的点可能被选作支配点)
			{//贪心地选择最佳节点作为新的支配点
					maxcoverage=Coverage[i];
					bestvertex=i;
			}
        if(maxcoverage==0)   //表示新添加任何点作为支配点都不会增加任何未满足点的满足度. 可以证明在无孤立点的图中这种情况不会发生
		{
			printf("No solution to the instance!\n");
			return 0;
		}
		DSList[bestvertex]=1;          //加入新的支配点
		dssize++;                      //支配点个数更新
		for(j=1;j<=gvertexnum;j++)     //更新该加入正影响集的点的各邻点的未满意度，只需要更新已经加入正影响集的领点就ok。
			if(GAdjMatrix[bestvertex][j])
			{    Satisfied[j]++;
			     if(Satisfied[j]==0)   //如果点j被支配次数现在能达到其点度一半,则归于已满足的点类，顺带如果他的邻点属于满意点。则该点加入满意点的集合。
				 {
					 uncoverednum--;
					 for(k=1;k<=gvertexnum;k++)                   //更新相关点的覆盖面Coverage[.]
							if(GAdjMatrix[j][k] && DSList[k]==0)  //邻点是非支配点时才有必要更新
								Coverage[k]--;
				 }
			}
	}


	//***用于调试阶段: 检查所得正影响的支配集是否正确
	if(!CheckPIDS(dssize))
	{	printf("There is some wrong with Greedy1_PIDS!\n");
		getch();
	}//***用于调试阶段

    printf("\nGreedy1_PIDS所求正影响支配集大小为(优化前): %5d\n",dssize);
    //优化阶段: 由所得正影响支配集导出一个极小正影响支配集
    int newdssize;
	newdssize=RefinePIDS(dssize);

	//***用于调试阶段: 检查所得正影响的支配集是否正确
	if(!CheckPIDS(newdssize))
	{	printf("There is some wrong with Greedy1_PIDS!\n");
		getch();
	}//***用于调试阶段

	printf("\nGreedy2_PIDS所求正影响支配集大小为(优化后): %5d\n",newdssize);
	return newdssize;//返回所求支配集大小
}




//***用于调试阶段: 检查一个正影响支配集的正确性***
bool CheckPIDS(int dssize)
{//用于调试阶段: 对任何算法得到的正影响支配集DSList[1...gvertexnum](其大小为dssize)进行正确性检查, 正确返回1; 否则返回0.
    int i,j;
	int newdssize;
	int temp;

	//检查每点邻域内是否有至少一半的支配点, 并统计实际支配点的个数newdssize
	newdssize=0;
	for(i=1;i<=gvertexnum;i++)
	{   if(DSList[i]==1)
	        newdssize++;
	    temp=0;
		for(j=1;j<=gvertexnum;j++)
			if(GAdjMatrix[i][j] && DSList[j]==1)
				temp++;
		if(temp<(DegreeList[i]+1)/2)  //如果任何点i的邻域内支配点个数不到其点度数的一半, 表示有错
			return 0;
	}
	if(newdssize!=dssize)             //如果实际的支配点个数与原来不符合, 表示有错
		return 0;
	return 1;
}





//****由一个正影响支配集导出一个极小正影响支配集****
int RefinePIDS(int dssize)
{//可对任何算法得到的正影响支配集DSList[1...gvertexnum]进行极小化处理(其大小为dssize),去掉不必要的支配点最终得到更小的正影响支配集DSList[1...gvertexnum],并返回新的大小
    int i,j,k;
	int newdssize;
	newdssize=dssize;  //获取原正影响支配集的大小

	//按照节点度由小到大的次序来进行下述操作(注意: 图中节点已经按照度由大到小有序,即第1个节点度最大,第gvertexnum节点度最小)
	for(i=gvertexnum;i>=1;i--)
		if(DSList[i]==1)  //表示DSList[i]是正影响集中一员，以此检查
		{	for(j=1;j<=gvertexnum;j++)
				if(GAdjMatrix[i][j] && Satisfied[j]==0) //此时说明支配点i不可以改为非支配点
					break;
			if(j==gvertexnum+1)                           //此时说明支配点i的任何邻点j都有Satisfied[j]<0, 因此节点i可以改为非支配点
			{
				newdssize--;
				DSList[i]=0;
		        for(k=1;k<=gvertexnum;k++)
			        if(GAdjMatrix[i][k])
			           Satisfied[k]--;
			}
		}

	return newdssize;//返回极小化后的正影响支配集的大小
}





/*现在开始复现师兄论文。
第一步，实现Crossver operator*/



void  CrossverOperator()
{

/*问题，这个x，y如何定制？*/

    bool      DSList_Z[GMaxVertexNum+1];       //存放Z
    bool      DSList_X[GMaxVertexNum+1];      //存放X。
    bool      DSList_Y[GMaxVertexNum+1];      //存放Y。
//如果随机生成x，y,一半是x，一半是y。
    for(int i=0; i<(strlen(DSList)/2); i++)
    {
        DSList_X[i]=DSList[i];
        DSList_Y[i]=DSList[(strlen(DSList)/2)+i];
    }

    //将T置为空集。
    for(int i=0;i<(strlen(DSList);i++){
         DSList[i]=0;
    }
    /*将x与y进行交叉。生成z。并入T*/
    for(int i=0; i<strlen(DSList_X); i++)
    {

        if(DSList_X[i]==DSList_Y[i]&&DSList_X[i]==1)
        {
            DSList_Z[i]=1;

        }
        DSList_Z[i]=0;
    }

    /*对z进行操作。*/

    for(int i=0; i<strlen(DSList_X); i++)
    {
        if(DSList_Z[i]==1)
        {
            //产生一个0到1的随机数。
            float RangNum=rand()/(RAND_MAX+1.0);
            if(RangNum<0.15)
            {
                DSList_Z[i]=0;

            }
        }
    }


   printf("\n产生的Z的数量: %5d\n",strlen(DSList_Z));
}





































//***文献中已有的贪心算法2***
int Greedy2_PIDS(void)
{//策略: 优先添加能提升其所有邻点未满足程度和最大(即被需要程度最大)的一个非支配点变为支配点,直到所有点都满足为止.
 //实现方法: 优先选取数组Needed[...]中取值最小(即绝对值最大)的那个非支配节点(DSList[.]值0的点是非支配点)变为支配点.
	int i,j,k,dssize;

	for(i=1;i<=gvertexnum;i++)
	{//初始化
		DSList[i]=0;                       //***记录支配点(值为0的点为非支配点,值为1的点为支配点)
		Satisfied[i]=-(DegreeList[i]+1)/2; //***记录每个点的未满足程度(值>=0表示是已满足点; 值<0表示未满足点,且其绝对值表示其未满足的程度)
	}

    for(i=1;i<=gvertexnum;i++)
	{//初始化
		Needed[i]=0;         //***记录每个非支配点(即DSList[.]值为0的点)的所有未满足邻点的未满足度之和(为负数,其绝对值表示该非支配点被需要的程度)
		for(j=1;j<=gvertexnum;j++)
           if(GAdjMatrix[i][j])
			   Needed[i]+=Satisfied[j];

	}
    dssize=0;//支配点个数初始化


	//预处理(略): 所有点度1的点邻点必定选作支配点(某些度为1可能也要被选作支配点,比如一个点与多个度1的点相邻时)
	//贪心算法
	int minneeded, bestvertex;
	int uncoverednum;
	uncoverednum=gvertexnum;   //初始化,表示开始时未满足点的个数

	while(uncoverednum)        //如果还有未满足的节点就需要继续添加新的支配点
	{
		minneeded=0;           //初始化为0.一般为负数
		for(i=1;i<=gvertexnum;i++)
			if(DSList[i]==0 && Needed[i]<minneeded)  //所有非支配点都可能要选取(包括某些度为1的点可能被选作支配点)
			{//贪心地选择最佳节点作为新的支配点
					minneeded=Needed[i];
					bestvertex=i;
			}
        if(minneeded==0)   //表示新添加任何点作为支配点都不会增加任何未满足点的满足度. 可以证明在无孤立点的图中这种情况不会发生
		{
			printf("No solution to the instance!\n");
			return 0;
		}
		DSList[bestvertex]=1;          //加入新的支配点
		dssize++;                      //支配点个数更新
		for(j=1;j<=gvertexnum;j++)     //更新各邻点的未满足程度
			if(GAdjMatrix[bestvertex][j])
			{    Satisfied[j]++;
  		         if(Satisfied[j]<=0)   //只有对之前的未满足的邻点(刚刚获得满足的点之前也是未满足点)才有必要更新相关信息
				 {	 for(k=1;k<=gvertexnum;k++)               //更新相关点的邻点的不满足之和Needed[]
    					if(GAdjMatrix[j][k] && DSList[k]==0)  //邻点是非支配点时才有必要更新
								Needed[k]++;                  //注: 值为负,其绝对值表示被需要的程度
				 }
				 if(Satisfied[j]==0)                          //如果点j被支配次数现在能达到其点度一半,则归于已满足的点类
						 uncoverednum--;
			}
	}

    //***用于调试阶段: 检查所得正影响的支配集是否正确
	if(!CheckPIDS(dssize))
	{	printf("There is some wrong with Greedy2_PIDS!\n");
		getch();
	}//***用于调试阶段

	printf("\nGreedy2_PIDS所求正影响支配集大小为(优化前): %5d\n",dssize);
    //优化阶段: 由所得正影响支配集导出一个极小正影响支配集
    int newdssize;
	newdssize=RefinePIDS(dssize);

	//***用于调试阶段: 检查所得正影响的支配集是否正确
	if(!CheckPIDS(newdssize))
	{	printf("There is some wrong with Greedy2_PIDS!\n");
		getch();
	}//***用于调试阶段

	printf("\nGreedy2_PIDS所求正影响支配集大小为(优化后): %5d\n",newdssize);
	return newdssize;//返回所求支配集大小
}


//***新设计的贪心算法1***
int NewGreedy1_PIDS(void)
{//该算法实际上是文献中两个贪心算法策略的一种结合的应用.
 //策略1:优先选择能提升最多个未满足点的满足度(即覆盖面最大)的一个非支配点变为支配点.
 //实现方法: 优先选取Coverage[...]中取最大值的那个非支配节点u(DSList[u]=0)变为支配点.
 //策略2: 如果有多个覆盖面最大的节点u的时, 优先选择那个能提升其邻点不满足程度和绝对值最大(即被需要程度最大)的一个非支配点u变为支配点.
 //实现方法: 在Coverage[u]值最大的多个节点u中优先选取Needed[u]值最小(即绝对值最大)的那个非支配节点u(DSList[u]=0)变为支配点.
	int i,j,k,dssize;

	for(i=1;i<=gvertexnum;i++)
	{//初始化
		DSList[i]=0;                         //***记录支配点(值为0的点为非支配点,值为1的点为支配点)
		Satisfied[i]=-(DegreeList[i]+1)/2;   //***记录每个点的未满足程度(值>=0表示是已满足点; 值<0表示未满足点,且其绝对值表示其未满足的程度)
		Coverage[i]=DegreeList[i];           //***记录每个非支配点(即DSList[.]值为0的点)的邻域内包含有的未满足的点个数--即记录每个非支配点的覆盖面
	}

    for(i=1;i<=gvertexnum;i++)
	{//初始化
		Needed[i]=0;                 //***记录每个非支配点(即DSList[.]值为0的点)的所有未满足邻点的未满足度之和(为负数,其绝对值表示该非支配点被需要的程度)
		for(j=1;j<=gvertexnum;j++)
           if(GAdjMatrix[i][j])
			   Needed[i]+=Satisfied[j];

	}
    dssize=0;//支配点个数初始化


	//预处理(略): 所有点度1的点邻点必定选作支配点(某些度为1可能也要被选作支配点,比如一个点与多个度1的点相邻时)
	//改进的贪心近似算法
	int maxcoverage, minneeded, bestvertex;
	int uncoverednum;
	uncoverednum=gvertexnum;   //初始化,表示开始时没有被满足点的个数

	while(uncoverednum)        //如果还有没有被满足的节点就需要继续添加新的支配点
	{
		maxcoverage=0;         //初始化.用于贪心策略1
		minneeded=0;           //初始化.用于贪心策略2
		for(i=1;i<=gvertexnum;i++)
		{	if(DSList[i]==0)   //所有非支配点都可能要选取(包括某些度为1的点可能须选作支配点)
			{
				if(Coverage[i]>maxcoverage)  //贪心策略1: 优先选择覆盖面大的非支配点为新的支配点
				{
					maxcoverage=Coverage[i];
					bestvertex=i;
				}
				else if(Coverage[i]==maxcoverage && Coverage[i]>0 && Needed[i]<minneeded)
				{//贪心策略2: 覆盖面同样最大的情况下则优先选择被需要程度最大(即Satisfied[u]最小的u)的非支配点为新的支配点
					minneeded=Needed[i];
					bestvertex=i;
				}
			}
		}
        if(maxcoverage==0)   //表示新添加任何点作为支配点都不会增加任何没有满足的点的被支配程度. 可以证明这种情况应该不会发生
		{
			printf("No solution to the instance!\n");
			return 0;
		}
		DSList[bestvertex]=1;          //加入新的支配点
		dssize++;                      //支配点个数更新
		for(j=1;j<=gvertexnum;j++)     //更新各邻点未满足程度
			if(GAdjMatrix[bestvertex][j])
			{    Satisfied[j]++;
				 if(Satisfied[j]<=0) //只有对之前的未满足的邻点(刚刚获得满足的点之前也是未满足点)才有必要更新相关信息
				 {	 for(k=1;k<=gvertexnum;k++)               //更新相关点的邻点的不满足之和Needed[.]**************
    					if(GAdjMatrix[j][k] && DSList[k]==0)  //邻点是非支配点时才有必要更新*************
								Needed[k]++;          //注意：值为负，其绝对值表示被需要的程度*************
				 }
			     if(Satisfied[j]==0) //如果点j被支配次数现在能达到其点度一半,则归于满足的点类
				 {	 uncoverednum--;
					 for(k=1;k<=gvertexnum;k++)               //更新相关点的覆盖面Coverage[.]**************
						if(GAdjMatrix[j][k] && DSList[k]==0)  //*************
							Coverage[k]--;                    //*************
				 }
			}
	}

   //***用于调试阶段: 检查所得正影响的支配集是否正确
	if(!CheckPIDS(dssize))
	{	printf("There is some wrong with NewGreedy1_PIDS!\n");
		getch();
	}//***用于调试阶段


    printf("\nNewGreedy1_PIDS所求正影响支配集大小为(优化前): %5d\n",dssize);
	//优化阶段: 由所得正影响支配集导出一个极小正影响支配集
    int newdssize;
	newdssize=RefinePIDS(dssize);

	//***用于调试阶段: 检查所得正影响的支配集是否正确
	if(!CheckPIDS(newdssize))
	{	printf("There is some wrong with NewGreedy1_PIDS!\n");
		getch();
	}//***用于调试阶段

	printf("\nNewGreedy1_PIDS所求正影响支配集大小为(优化后): %5d\n",newdssize);
	return newdssize;//返回所求支配集大小
}



//***新设计的贪心算法2***
int NewGreedy2_PIDS(void)
{//该算法实际上是文献中两个贪心算法策略的一种结合的应用.
 //策略2: 优先选择能提升其邻点不满足程度和绝对值最大(即被需要程度最大)的一个非支配点u变为支配点.
 //实现方法: 优先选取Needed[u]值最小(即绝对值最大)的那个非支配节点u(DSList[u]=0)变为支配点.
 //策略1: 如果有多个被需要程度最大的节点时, 则优先选择能提升最多个未满足点的满足度(即覆盖面最大)的一个非支配点u变为支配点.
 //实现方法: 在Needed[u]绝对值最大的多个节点u中优先选取Coverage[...]最大值的那个非支配节点u(DSList[u]=0)变为支配点.

	int i,j,k,dssize;

	for(i=1;i<=gvertexnum;i++)
	{//初始化
		DSList[i]=0;                         //***记录支配点(值为0的点为非支配点,值为1的点为支配点)
		Satisfied[i]=-(DegreeList[i]+1)/2;   //***记录每个点的未满足程度(值>=0表示是已满足点; 值<0表示未满足点,且其绝对值表示其未满足的程度)
		Coverage[i]=DegreeList[i];           //***记录每个非支配点(即DSList[.]值为0的点)的邻域内包含有的未满足的点个数--即记录每个非支配点的覆盖面
	}

    for(i=1;i<=gvertexnum;i++)
	{//初始化
		Needed[i]=0;                 //***记录每个非支配点(即DSList[.]值为0的点)的所有未满足邻点的未满足度之和(为负数,其绝对值表示该非支配点被需要的程度)
		for(j=1;j<=gvertexnum;j++)
           if(GAdjMatrix[i][j])
			   Needed[i]+=Satisfied[j];

	}
    dssize=0;//支配点个数初始化


	//预处理(略): 所有点度1的点邻点必定选作支配点(某些度为1可能也要被选作支配点,比如一个点与多个度1的点相邻时)
	//改进的贪心近似算法
	int maxcoverage, minneeded, bestvertex;
	int uncoverednum;
	uncoverednum=gvertexnum;   //初始化,表示开始时没有被满足点的个数

	while(uncoverednum)        //如果还有没有被满足的节点就需要继续添加新的支配点
	{
		maxcoverage=0;         //初始化.用于贪心策略1
		minneeded=0;           //初始化.用于贪心策略2
		for(i=1;i<=gvertexnum;i++)
		{	if(DSList[i]==0)   //所有非支配点都可能要选取(包括某些度为1的点可能须选作支配点)
			{
				if(minneeded>Needed[i])  //贪心策略2: 优先选择被需要程度最大(即Needed[u]最小的u)的非支配点为新的支配点
				{
					minneeded=Needed[i];
					bestvertex=i;
				}
				else if(Needed[i]==minneeded && Needed[i]<0 && Coverage[i]>maxcoverage)
				{//贪心策略1: 被需要程度同样最大的情况下则优先选择覆盖面大的一个非支配点为新的支配点
					maxcoverage=Coverage[i];
					bestvertex=i;
				}
			}
		}
        if(minneeded==0)   //表示新添加任何点作为支配点都不会增加任何没有满足的点的被支配程度. 可以证明这种情况应该不会发生
		{
			printf("No solution to the instance!\n");
			return 0;
		}
		DSList[bestvertex]=1;          //加入新的支配点
		dssize++;                      //支配点个数更新
		for(j=1;j<=gvertexnum;j++)     //更新各邻点未满足程度
			if(GAdjMatrix[bestvertex][j])
			{    Satisfied[j]++;
				 if(Satisfied[j]<=0)   //只有对之前的未满足的邻点(刚刚获得满足的点之前也是未满足点)才有必要更新相关信息
				 {	 for(k=1;k<=gvertexnum;k++)               //更新相关点的邻点的不满足之和Needed[.]**************
    					if(GAdjMatrix[j][k] && DSList[k]==0)  //邻点是非支配点时才有必要更新*************
								Needed[k]++;                  //注：值为负，其绝对值表示被需要的程度*************
				 }
			     if(Satisfied[j]==0) //如果点j被支配次数现在能达到其点度一半,则归于满足的点类
				 {	 uncoverednum--;
					 for(k=1;k<=gvertexnum;k++)               //更新相关点的覆盖面Coverage[.]**************
						if(GAdjMatrix[j][k] && DSList[k]==0)  //*************
							Coverage[k]--;                    //*************
				 }
			}
	}

   //***用于调试阶段: 检查所得正影响的支配集是否正确
	if(!CheckPIDS(dssize))
	{	printf("There is some wrong with NewGreedy2_PIDS!\n");
		getch();
	}//***用于调试阶段

    printf("\nNewGreedy2_PIDS所求正影响支配集大小为(优化前): %5d\n",dssize);
    //优化阶段: 由所得正影响支配集导出一个极小正影响支配集
    int newdssize;
	newdssize=RefinePIDS(dssize);

	//***用于调试阶段: 检查所得正影响的支配集是否正确
	if(!CheckPIDS(newdssize))
	{	printf("There is some wrong with NewGreedy2_PIDS!\n");
		getch();
	}//***用于调试阶段

	printf("\nNewGreedy2_PIDS所求正影响支配集大小为(优化后): %5d\n",newdssize);
	return newdssize;//返回所求支配集大小
}


//***新设计的贪心算法3***
int NewGreedy3_PIDS(void)
{//该算法实际上是文献中两个贪心算法策略的一种结合的应用.
 //策略2: 优先选择能提升其邻点不满足程度和绝对值最大(即被需要程度最大)的一个非支配点u变为支配点.
 //实现方法: 优先选取Needed[u]值最小(即绝对值最大)的那个非支配节点u(DSList[u]=0)变为支配点.
 //策略1: 如果有多个被需要程度最大的节点时, 则优先选择能提升最小个未满足点的满足度(即覆盖面最 小)的一个非支配点u变为支配点.
 //实现方法: 在Needed[u]绝对值最大的多个节点u中优先选取Coverage[...]最小 值(此时平均每个未满足点的不满足度很大)的那个非支配节点u(DSList[u]=0)变为支配点.

	int i,j,k,dssize;

	for(i=1;i<=gvertexnum;i++)
	{//初始化
		DSList[i]=0;                         //***记录支配点(值为0的点为非支配点,值为1的点为支配点)
		Satisfied[i]=-(DegreeList[i]+1)/2;   //***记录每个点的未满足程度(值>=0表示是已满足点; 值<0表示未满足点,且其绝对值表示其未满足的程度)
		Coverage[i]=DegreeList[i];           //***记录每个非支配点(即DSList[.]值为0的点)的邻域内包含有的未满足的点个数--即记录每个非支配点的覆盖面
	}

    for(i=1;i<=gvertexnum;i++)
	{//初始化
		Needed[i]=0;                 //***记录每个非支配点(即DSList[.]值为0的点)的所有未满足邻点的未满足度之和(为负数,其绝对值表示该非支配点被需要的程度)
		for(j=1;j<=gvertexnum;j++)
           if(GAdjMatrix[i][j])
			   Needed[i]+=Satisfied[j];

	}
    dssize=0;//支配点个数初始化


	//预处理(略): 所有点度1的点邻点必定选作支配点(某些度为1可能也要被选作支配点,比如一个点与多个度1的点相邻时)
	//改进的贪心近似算法
	int mincoverage, minneeded, bestvertex;
	int uncoverednum;
	uncoverednum=gvertexnum;   //初始化,表示开始时没有被满足点的个数

	while(uncoverednum)        //如果还有没有被满足的节点就需要继续添加新的支配点
	{
		mincoverage=gvertexnum;  //初始化.用于贪心策略1
		minneeded=0;             //初始化.用于贪心策略2
		for(i=1;i<=gvertexnum;i++)
		{	if(DSList[i]==0)   //所有非支配点都可能要选取(包括某些度为1的点可能须选作支配点)
			{
				if(minneeded>Needed[i])  //贪心策略2: 优先选择被需要程度最大(即Needed[u]最小的u)的非支配点为新的支配点
				{
					minneeded=Needed[i];
					bestvertex=i;
				}
				else if(Needed[i]==minneeded && Needed[i]<0 && Coverage[i]<mincoverage)
				{//贪心策略1: 被需要程度同样最大的情况下则优先选择覆盖面大的一个非支配点为新的支配点
					mincoverage=Coverage[i];
					bestvertex=i;
				}
			}
		}
        if(minneeded==0)   //表示新添加任何点作为支配点都不会增加任何没有满足的点的被支配程度. 可以证明这种情况应该不会发生
		{
			printf("No solution to the instance!\n");
			return 0;
		}
		DSList[bestvertex]=1;          //加入新的支配点
		dssize++;                      //支配点个数更新
		for(j=1;j<=gvertexnum;j++)     //更新各邻点未满足程度
			if(GAdjMatrix[bestvertex][j])
			{    Satisfied[j]++;
				 if(Satisfied[j]<=0)   //只有对之前的未满足的邻点(刚刚获得满足的点之前也是未满足点)才有必要更新相关信息
				 {	 for(k=1;k<=gvertexnum;k++)               //更新相关点的邻点的不满足之和Needed[.]**************
    					if(GAdjMatrix[j][k] && DSList[k]==0)  //邻点是非支配点时才有必要更新*************
								Needed[k]++;                  //注：值为负，其绝对值表示被需要的程度*************
				 }
			     if(Satisfied[j]==0) //如果点j被支配次数现在能达到其点度一半,则归于满足的点类
				 {	 uncoverednum--;
					 for(k=1;k<=gvertexnum;k++)               //更新相关点的覆盖面Coverage[.]**************
						if(GAdjMatrix[j][k] && DSList[k]==0)  //*************
							Coverage[k]--;                    //*************
				 }
			}
	}

   //***用于调试阶段: 检查所得正影响的支配集是否正确
	if(!CheckPIDS(dssize))
	{	printf("There is some wrong with NewGreedy3_PIDS!\n");
		getch();
	}//***用于调试阶段

    printf("\nNewGreedy3_PIDS所求正影响支配集大小为(优化前): %5d\n",dssize);
    //优化阶段: 由所得正影响支配集导出一个极小正影响支配集
    int newdssize;
	newdssize=RefinePIDS(dssize);

	//***用于调试阶段: 检查所得正影响的支配集是否正确
	if(!CheckPIDS(newdssize))
	{	printf("There is some wrong with NewGreedy3_PIDS!\n");
		getch();
	}//***用于调试阶段

	printf("\nNewGreedy3_PIDS所求正影响支配集大小为(优化后): %5d\n",newdssize);
	return newdssize;//返回所求支配集大小
}



//***新设计的局部贪心方法1***
int LocalGreedy1_PIDS(void)
{//策略: 每次选择一个最不满足点u, 然后在u的非支配点型邻点中反复使用策略1直到u满足为止.
 //最不满足点u是Satisfied[.]值为负且绝对值最大的点.在点u的邻域中优先依次选择覆盖面最大的若干非支配点变为支配点(即优先选择Coverage[.]值大且与u相邻的非支配点),直到点u满足为止.

	int i,j,k,dssize;

	for(i=1;i<=gvertexnum;i++)
	{//初始化
		DSList[i]=0;                       //***记录支配点(值为0的点为非支配点,值为1的点为支配点)
		Satisfied[i]=-(DegreeList[i]+1)/2; //***记录每个点的未满足程度(值>=0表示是已满足点; 值<0表示未满足点,且其绝对值表示其未满足的程度)
		Coverage[i]=DegreeList[i];         //***记录每个非支配点(即DSList[.]值为0的点)的邻域内包含有的未满足的点个数--即记录每个非支配点的覆盖面
	}
    dssize=0;//支配点个数初始化

	//局部贪心算法
	int uncoverednum;
	int currentvertex, minsatisfied, maxcoverage, maxcvertex;

	uncoverednum=gvertexnum;   //初始化,表示开始时未满足点的个数
	while(uncoverednum)
	{//如果uncoverednum>0, 即还有未满足的节点就需要继续添加新的支配点

        minsatisfied=0;        //用于找出满足程度最小的节点(即Satisfied[]值为负且值最小的点,其绝对值表示其未满足程度的最大)
		currentvertex=0;
		for(i=1;i<=gvertexnum;i++)   //寻找当前满足程度最小的节点处理(如果简单点处理,也可以按点度最大的未满足节点优先代替)
		    if(Satisfied[i]<minsatisfied)
			{
		    	minsatisfied=Satisfied[i];
		    	currentvertex=i;
			}

		if(currentvertex==0) break;

        //******局部贪心策略:一次性将当前未满足的节点currentvertex变为满足节点(依次将其邻域中覆盖面最大的若干非支配点变为支配点)
		while(Satisfied[currentvertex]<0) //Satisfied[currentvertex]<0表示currentvertex是未满足点
		{
			maxcoverage=0;   //用于找出能覆盖未满足节点最多的一个非支配点(currentvertex的邻点)

			//可以取消下面几行标记为  //**的代码(如果实验发现不考虑当前点也可变为支配点时效果更好的话)
			//if(DSList[currentvertex]==0 && Coverage[currentvertex]>maxcoverage)  //** 如果当前未满足节点currentvertex是非支配点时,也考虑它可变为支配点
			//{	maxcoverage=Coverage[currentvertex];  //*
			//	maxcvertex=currentvertex;      //*
			//}

			for(i=1;i<=gvertexnum;i++)
			{
				if(GAdjMatrix[currentvertex][i] && DSList[i]==0 && Coverage[i]>maxcoverage)  //局部贪心地选择一个邻点(非支配点)
				{   	maxcoverage=Coverage[i];
						maxcvertex=i;
				}
			}
		    if(maxcoverage==0)   //表示新添加任何邻点作为支配点都不会增加任何未满足的点的满足程度. 可以证明这种情况在不含孤立点的图中不会发生
			{
				printf("No solution to the instance!\n");
				return 0;
			}

			DSList[maxcvertex]=1;          //加入新的支配点
			dssize++;                      //支配点个数更新
			for(j=1;j<=gvertexnum;j++)     //更新相关点的未满足程度Satisfied[.]
				if(GAdjMatrix[maxcvertex][j])
				{   Satisfied[j]++;
					if(Satisfied[j]==0)    //如果点j被支配次数现在能达到其点度一半,则归于满足的点类
					{	uncoverednum--;
					    for(k=1;k<=gvertexnum;k++)               //更新相关点的覆盖面Coverage[.]**************
							if(GAdjMatrix[j][k] && DSList[k]==0)
								Coverage[k]--;
					}
				}
		}
	}


	//***用于调试阶段: 检查所得正影响的支配集是否正确
	if(!CheckPIDS(dssize))
	{	printf("There is some wrong with LocalGreedy1_PIDS!\n");
		getch();
	}//***用于调试阶段

	printf("\nLocalGreedy1_PIDS所求正影响支配集大小为(优化前): %5d\n",dssize);
	//优化阶段: 由所得正影响支配集导出一个极小正影响支配集
    int newdssize;
	newdssize=RefinePIDS(dssize);

	//***用于调试阶段: 检查所得正影响的支配集是否正确
	if(!CheckPIDS(newdssize))
	{	printf("There is some wrong with LocaGreedy1_PIDS!\n");
		getch();
	}//***用于调试阶段

	printf("\nLocalGreedy1_PIDS所求正影响支配集大小为(优化后): %5d\n",newdssize);
	return newdssize;//返回所求支配集大小
}



//***新设计的局部贪心方法2***
int LocalGreedy2_PIDS(void)
{//策略: 每次选择一个最不满足点u, 然后在u的非支配点型邻点中反复使用策略2直到u满足为止.
 //最不满足点u是Satisfied[.]值为负且绝对值最大的点.在点u的邻域中优先依次选择被需要程度最大的若干非支配点变为支配点(即优先选择Needed[.]绝对值大且与u相邻的非支配点),直到点u满足为止.

	int i,j,k,dssize;

	for(i=1;i<=gvertexnum;i++)
	{//初始化
		DSList[i]=0;                       //***记录支配点(值为0的点为非支配点,值为1的点为支配点)
		Satisfied[i]=-(DegreeList[i]+1)/2; //***记录每个点的未满足程度(值>=0表示是已满足点; 值<0表示未满足点,且其绝对值表示其未满足的程度)
	}
    for(i=1;i<=gvertexnum;i++)
	{//初始化
		Needed[i]=0;         //***记录每个非支配点(即DSList[.]值为0的点)的所有未满足邻点的未满足度之和(为负数,其绝对值表示该非支配点被需要的程度)
		for(j=1;j<=gvertexnum;j++)
           if(GAdjMatrix[i][j])
			   Needed[i]+=Satisfied[j];

	}
    dssize=0;//支配点个数初始化

	//局部贪心算法
	int uncoverednum;
	int currentvertex, minsatisfied, minneeded, minneededvertex;

	uncoverednum=gvertexnum;         //初始化,表示开始时未满足点的个数
	while(uncoverednum)
	{//如果uncoverednum>0, 即还有未满足的节点就需要继续添加新的支配点

        minsatisfied=0;              //用于找出满足程度最小的节点(即Satisfied[]值为负且值最小的点,其绝对值表示其未满足程度的最大)
		currentvertex=0;
		for(i=1;i<=gvertexnum;i++)   //寻找当前满足程度最小的节点处理(如果简单点处理,也可以按点度最大的未满足节点优先代替)
		    if(Satisfied[i]<minsatisfied)
			{
		    	minsatisfied=Satisfied[i];
		    	currentvertex=i;
			}

		if(currentvertex==0) break;

        //******局部贪心法策略:一次性将当前未满足的节点currentvertex变为满足节点(依次将其邻域中被需要程度最大的若干非支配点变为支配点)
		while(Satisfied[currentvertex]<0) //Satisfied[currentvertex]<0表示currentvertex是未满足点
		{
			minneeded=0;   //用于找出被需要程度最大(即Needed[]值为负且绝对值最大)的一个非支配点(currentvertex的邻点)

			//可以取消下面几行标记为  //**的代码(如果实验发现不考虑当前点也可变为支配点时效果更好的话)
			//if(DSList[currentvertex]==0 && Needed[currentvertex]<minneeded)  //** 如果当前未满足节点currentvertex是非支配点时,也考虑它可变为支配点
			//{	minneeded=Needed[currentvertex];  //*
			//	minneededvertex=currentvertex;      //*
			//}

			for(i=1;i<=gvertexnum;i++)
			{
				if(GAdjMatrix[currentvertex][i] && DSList[i]==0 && Needed[i]<minneeded)  //局部贪心地选择一个邻点(非支配点)
				{   	minneeded=Needed[i];
						minneededvertex=i;
				}
			}
		    if(minneeded==0)   //表示新添加任何邻点作为支配点都不会增加任何未满足的点的满足程度. 可以证明这种情况在不含孤立点的图中不会发生
			{
				printf("No solution to the instance!\n");
				return 0;
			}

			DSList[minneededvertex]=1;          //加入新的支配点
			dssize++;                      //支配点个数更新
			for(j=1;j<=gvertexnum;j++)     //更新相关点的未满足程度Satisfied[.]
				if(GAdjMatrix[minneededvertex][j])
				{    Satisfied[j]++;
  					if(Satisfied[j]<=0)   //只有对之前的未满足的邻点(刚刚获得满足的点之前也是未满足点)才有必要更新相关信息
					{	for(k=1;k<=gvertexnum;k++)                //更新相关点的邻点的不满足之和Needed[]
    						if(GAdjMatrix[j][k] && DSList[k]==0)  //邻点是非支配点时才有必要更新
								Needed[k]++;                      //注: 值为负,其绝对值表示被需要的程度
					}
					if(Satisfied[j]==0)                           //如果点j被支配次数现在能达到其点度一半,则归于已满足的点类
						 uncoverednum--;
				}
		}
	}

    //***用于调试阶段: 检查所得正影响的支配集是否正确
	if(!CheckPIDS(dssize))
	{	printf("There is some wrong with LocalGreedy2_PIDS!\n");
		getch();
	}//***用于调试阶段

	printf("\nLocalGreedy2_PIDS所求正影响支配集大小为(优化前): %5d\n",dssize);

    //优化阶段: 由所得正影响支配集导出一个极小正影响支配集
    int newdssize;
	newdssize=RefinePIDS(dssize);

	//***用于调试阶段: 检查所得正影响的支配集是否正确
	if(!CheckPIDS(newdssize))
	{	printf("There is some wrong with LocaGreedy2_PIDS!\n");
		getch();
	}//***用于调试阶段

	printf("\nLocalGreedy2_PIDS所求正影响支配集大小为(优化后): %5d\n",newdssize);
	return newdssize;//返回所求支配集大小
}


//***新设计的局部贪心方法3(快速算法--该算法时间复杂度为O(n^2),其他算法时间复杂都是O(n^3))***
int LocalGreedy3_PIDS(void)
{//策略: 每次选择一个最不满足点u, 然后在u的非支配点型邻点中反复添加点度大的点作为新的支配点,直到u满足为止.
 //最不满足点u是Satisfied[.]值为负且绝对值最大的点.点u邻域中依次添加点度大的点作为新的支配点(读图或建图时预处理中图点序已按点度降序排了,故自然次序即可),直到点u满足为止.

	int i,j,dssize;

	for(i=1;i<=gvertexnum;i++)
	{//初始化
		DSList[i]=0;                       //***记录支配点(值为0的点为非支配点,值为1的点为支配点)
		Satisfied[i]=-(DegreeList[i]+1)/2; //***记录每个点的未满足程度(值>=0表示是已满足点; 值<0表示未满足点,且其绝对值表示其未满足的程度)
	}

    dssize=0;//支配点个数初始化

	//局部贪心算法
	int uncoverednum;
	int currentvertex, minsatisfied, bestvertex;

	uncoverednum=gvertexnum;         //初始化,表示开始时未满足点的个数
	while(uncoverednum)
	{//如果uncoverednum>0, 即还有未满足的节点就需要继续添加新的支配点

        minsatisfied=0;              //用于找出满足程度最小的节点(即Satisfied[]值为负且值最小的点,其绝对值表示其未满足程度的最大)
		currentvertex=0;
		for(i=1;i<=gvertexnum;i++)   //寻找当前满足程度最小的节点处理(如果简单点处理,也可以按点度最大的未满足节点优先代替)
		    if(Satisfied[i]<minsatisfied)
			{
		    	minsatisfied=Satisfied[i];
		    	currentvertex=i;
			}

		if(currentvertex==0) break;

        //******局部贪心策略:一次性将当前未满足的节点currentvertex变为满足节点(依次将其邻域中点度最大的若干非支配点变为支配点)
		while(Satisfied[currentvertex]<0) //Satisfied[currentvertex]<0表示currentvertex是未满足点
		{
			for(i=1;i<=gvertexnum;i++)
				if(GAdjMatrix[currentvertex][i] && DSList[i]==0)  //局部贪心地选择一个邻点(非支配点)--点度大优先(自然次序即表示度大的优先,这是因为图点序已按点度降序排了)
				{		bestvertex=i;
				        break;
				}
		    if(i==gvertexnum+1)   //表示新添加任何邻点作为支配点都不会增加任何未满足的点的满足程度. 可以证明这种情况在不含孤立点的图中不会发生
			{
				printf("No solution to the instance!\n");
				return 0;
			}

			DSList[bestvertex]=1;          //加入新的支配点
			dssize++;                      //支配点个数更新
			for(j=1;j<=gvertexnum;j++)     //更新相关点的未满足程度Satisfied[.]
				if(GAdjMatrix[bestvertex][j])
				{    Satisfied[j]++;
  					 if(Satisfied[j]==0)   //如果点j被支配次数现在能达到其点度一半,则归于已满足的点类
						 uncoverednum--;
				}
		}
	}


    //***用于调试阶段: 检查所得正影响的支配集是否正确
	if(!CheckPIDS(dssize))
	{	printf("There is some wrong with LocalGreedy3_PIDS!\n");
		getch();
	}//***用于调试阶段

	printf("\nLocalGreedy3_PIDS所求正影响支配集大小为(优化前): %5d\n",dssize);

    //优化阶段: 由所得正影响支配集导出一个极小正影响支配集
    int newdssize;
	newdssize=RefinePIDS(dssize);

	//***用于调试阶段: 检查所得正影响的支配集是否正确
	if(!CheckPIDS(newdssize))
	{	printf("There is some wrong with LocaGreedy3_PIDS!\n");
		getch();
	}//***用于调试阶段

	printf("\nLocalGreedy3_PIDS所求正影响支配集大小为(优化后): %5d\n",newdssize);
	return newdssize;//返回所求支配集大小
}




//***新设计的新的局部贪心方法1***
int NewLocalGreedy1_PIDS(void)
{//策略: 每次选择当前最不满足点u, 然后在u的非支配点型邻点中使用策略1(使用1次)
 //最不满足点u是Satisfied[.]值为负且绝对值最大的点.在点u的邻域中选择覆盖面最大(即Coverage[.]值大且与u相邻)的一个非支配点变为支配点.

	int i,j,k,dssize;

	for(i=1;i<=gvertexnum;i++)
	{//初始化
		DSList[i]=0;                       //***记录支配点(值为0的点为非支配点,值为1的点为支配点)
		Satisfied[i]=-(DegreeList[i]+1)/2; //***记录每个点的未满足程度(值>=0表示是已满足点; 值<0表示未满足点,且其绝对值表示其未满足的程度)
		Coverage[i]=DegreeList[i];         //***记录每个非支配点(即DSList[.]值为0的点)的邻域内包含有的未满足的点个数--即记录每个非支配点的覆盖面
	}
    dssize=0;//支配点个数初始化

	//新的局部贪心算法
	int uncoverednum;
	int currentvertex, minsatisfied, maxcoverage, maxcvertex;

	uncoverednum=gvertexnum;   //初始化,表示开始时未满足点的个数
	while(uncoverednum)
	{//如果uncoverednum>0, 即还有未满足的节点就需要继续添加新的支配点

        minsatisfied=0;        //用于找出满足程度最小的节点(即Satisfied[]值为负且值最小的点,其绝对值表示其未满足程度的最大)
		currentvertex=0;
		for(i=1;i<=gvertexnum;i++)   //寻找当前满足程度最小的节点处理(如果简单点处理,也可以按点度最大的未满足节点优先代替)
		    if(Satisfied[i]<minsatisfied)
			{
		    	minsatisfied=Satisfied[i];
		    	currentvertex=i;
			}

		if(currentvertex==0) break;

        //******新的局部贪心策略: 将当前未满足节点currentvertex的邻域中覆盖面最大的一个非支配点变为支配点
		maxcoverage=0;   //用于找出能覆盖未满足节点最多的一个非支配点(currentvertex的邻点)

		//可以取消下面几行标记为  //**的代码(如果实验发现不考虑当前点也可变为支配点时效果更好的话)
		//if(DSList[currentvertex]==0 && Coverage[currentvertex]>maxcoverage)  //** 如果当前未满足节点currentvertex是非支配点时,也考虑它可变为支配点
		//{	maxcoverage=Coverage[currentvertex];  //*
		//	maxcvertex=currentvertex;      //*
		//}

		for(i=1;i<=gvertexnum;i++)
			if(GAdjMatrix[currentvertex][i] && DSList[i]==0 && Coverage[i]>maxcoverage)  //局部贪心地选择一个邻点(非支配点)
			{   	maxcoverage=Coverage[i];
					maxcvertex=i;
			}
		if(maxcoverage==0)   //表示新添加任何邻点作为支配点都不会增加任何未满足的点的满足程度. 可以证明这种情况在不含孤立点的图中不会发生
		{
			printf("No solution to the instance!\n");
			return 0;
		}

		DSList[maxcvertex]=1;          //加入新的支配点
		dssize++;                      //支配点个数更新
		for(j=1;j<=gvertexnum;j++)     //更新相关点的未满足程度Satisfied[.]
			if(GAdjMatrix[maxcvertex][j])
			{   Satisfied[j]++;
				if(Satisfied[j]==0)    //如果点j被支配次数现在能达到其点度一半,则归于满足的点类
				{	uncoverednum--;
				    for(k=1;k<=gvertexnum;k++)               //更新相关点的覆盖面Coverage[.]**************
						if(GAdjMatrix[j][k] && DSList[k]==0)
							Coverage[k]--;
				}
			}
		}


	//***用于调试阶段: 检查所得正影响的支配集是否正确
	if(!CheckPIDS(dssize))
	{	printf("There is some wrong with NewLocalGreedy1_PIDS!\n");
		getch();
	}//***用于调试阶段

	printf("\nNewLocalGreedy1_PIDS所求正影响支配集大小为(优化前): %5d\n",dssize);
	//优化阶段: 由所得正影响支配集导出一个极小正影响支配集
    int newdssize;
	newdssize=RefinePIDS(dssize);

	//***用于调试阶段: 检查所得正影响的支配集是否正确
	if(!CheckPIDS(newdssize))
	{	printf("There is some wrong with NewLocaGreedy1_PIDS!\n");
		getch();
	}//***用于调试阶段

	printf("\nNewLocalGreedy1_PIDS所求正影响支配集大小为(优化后): %5d\n",newdssize);
	return newdssize;//返回所求支配集大小
}


//***新设计的新的局部贪心方法2***
int NewLocalGreedy2_PIDS(void)
{//策略: 每次选择当前最不满足点u, 然后在u的非支配点型邻点中使用策略2(使用1次).
 //最不满足点u是Satisfied[]值为负且绝对值最大的点.在点u邻域中选择被需要程度最大(即Needed[.]绝对值最大且与u相邻)的一个非支配点w变为支配点.
 //含义:u最需要别人且w也最被人需要
	int i,j,k,dssize;

	for(i=1;i<=gvertexnum;i++)
	{//初始化
		DSList[i]=0;                       //***记录支配点(值为0的点为非支配点,值为1的点为支配点)
		Satisfied[i]=-(DegreeList[i]+1)/2; //***记录每个点的未满足程度(值>=0表示是已满足点; 值<0表示未满足点,且其绝对值表示其未满足的程度)
	}
    for(i=1;i<=gvertexnum;i++)
	{//初始化
		Needed[i]=0;         //***记录每个非支配点(即DSList[.]值为0的点)的所有未满足邻点的未满足度之和(为负数,其绝对值表示该非支配点被需要的程度)
		for(j=1;j<=gvertexnum;j++)
           if(GAdjMatrix[i][j])
			   Needed[i]+=Satisfied[j];

	}
    dssize=0;//支配点个数初始化

	//局部贪心算法
	int uncoverednum;
	int currentvertex, minsatisfied, minneeded, minneededvertex;

	uncoverednum=gvertexnum;         //初始化,表示开始时未满足点的个数
	while(uncoverednum)
	{//如果uncoverednum>0, 即还有未满足的节点就需要继续添加新的支配点

        minsatisfied=0;              //用于找出满足程度最小的节点(即Satisfied[]值为负且值最小的点,其绝对值表示其未满足程度的最大)
		currentvertex=0;
		for(i=1;i<=gvertexnum;i++)   //寻找当前满足程度最小的节点处理(如果简单点处理,也可以按点度最大的未满足节点优先代替)
		    if(Satisfied[i]<minsatisfied)
			{
		    	minsatisfied=Satisfied[i];
		    	currentvertex=i;
			}

		if(currentvertex==0) break;

        //******局部贪心法策略:将当前未满足的节点currentvertex邻域中被需要程度最大的一个非支配点变为支配点
		minneeded=0;   //用于找出被需要程度最大(即Needed[]值为负且绝对值最大)的一个非支配点(currentvertex的邻点)

		//可以取消下面几行标记为  //**的代码(如果实验发现不考虑当前点也可变为支配点时效果更好的话)
		//if(DSList[currentvertex]==0 && Needed[currentvertex]<minneeded)  //** 如果当前未满足节点currentvertex是非支配点时,也考虑它可变为支配点
		//{	minneeded=Needed[currentvertex];  //*
		//	minneededvertex=currentvertex;      //*
		//}
		for(i=1;i<=gvertexnum;i++)
			if(GAdjMatrix[currentvertex][i] && DSList[i]==0 && Needed[i]<minneeded)  //局部贪心地选择一个邻点(非支配点)
			{   	minneeded=Needed[i];
					minneededvertex=i;
			}
	    if(minneeded==0)   //表示新添加任何邻点作为支配点都不会增加任何未满足的点的满足程度. 可以证明这种情况在不含孤立点的图中不会发生
		{
			printf("No solution to the instance!\n");
			return 0;
		}

		DSList[minneededvertex]=1;     //加入新的支配点
		dssize++;                      //支配点个数更新
		for(j=1;j<=gvertexnum;j++)     //更新相关点的未满足程度Satisfied[.]
			if(GAdjMatrix[minneededvertex][j])
			{    Satisfied[j]++;
 				 if(Satisfied[j]<=0)   //只有对之前的未满足的邻点(刚刚获得满足的点之前也是未满足点)才有必要更新相关信息
				 {	for(k=1;k<=gvertexnum;k++)                //更新相关点的邻点的不满足之和Needed[]
    					if(GAdjMatrix[j][k] && DSList[k]==0)  //邻点是非支配点时才有必要更新
							Needed[k]++;                      //注: 值为负,其绝对值表示被需要的程度
				}
				if(Satisfied[j]==0)                           //如果点j被支配次数现在能达到其点度一半,则归于已满足的点类
					 uncoverednum--;
			}
	}

    //***用于调试阶段: 检查所得正影响的支配集是否正确
	if(!CheckPIDS(dssize))
	{	printf("There is some wrong with NewLocalGreedy2_PIDS!\n");
		getch();
	}//***用于调试阶段

	printf("\nNewLocalGreedy2_PIDS所求正影响支配集大小为(优化前): %5d\n",dssize);

    //优化阶段: 由所得正影响支配集导出一个极小正影响支配集
    int newdssize;
	newdssize=RefinePIDS(dssize);

	//***用于调试阶段: 检查所得正影响的支配集是否正确
	if(!CheckPIDS(newdssize))
	{	printf("There is some wrong with NewLocaGreedy2_PIDS!\n");
		getch();
	}//***用于调试阶段

	printf("\nNewLocalGreedy2_PIDS所求正影响支配集大小为(优化后): %5d\n",newdssize);
	return newdssize;//返回所求支配集大小
}


//新设计的新的局部贪心方法3(快速算法--该算法时间复杂度为O(n^2),其它算法时间复杂都是O(n^3))
int NewLocalGreedy3_PIDS(void)
{//策略: 每次选择当前最不满足点u, 然后在u的非支配点型邻点中添加最度的点为新支配点.(读图或建图时预处理中图点序已按点度降序排了)
 //最不满足点u是Satisfied[.]值为负且绝对值最大的点.点u邻域中添加度最大的一个非支配型邻点作为新的支配点
	int i,j,dssize;

	for(i=1;i<=gvertexnum;i++)
	{//初始化
		DSList[i]=0;                       //***记录支配点(值为0的点为非支配点,值为1的点为支配点)
		Satisfied[i]=-(DegreeList[i]+1)/2; //***记录每个点的未满足程度(值>=0表示是已满足点; 值<0表示未满足点,且其绝对值表示其未满足的程度)
	}

    dssize=0;//支配点个数初始化

	//局部贪心算法
	int uncoverednum;
	int currentvertex, minsatisfied, bestvertex;

	uncoverednum=gvertexnum;         //初始化,表示开始时未满足点的个数
	while(uncoverednum)
	{//如果uncoverednum>0, 即还有未满足的节点就需要继续添加新的支配点

        minsatisfied=0;              //用于找出满足程度最小的节点(即Satisfied[]值为负且值最小的点,其绝对值表示其未满足程度的最大)
		currentvertex=0;
		for(i=1;i<=gvertexnum;i++)   //寻找当前满足程度最小的节点处理(如果简单点处理,也可以按点度最大的未满足节点优先代替)
		    if(Satisfied[i]<minsatisfied)
			{
		    	minsatisfied=Satisfied[i];
		    	currentvertex=i;
			}

		if(currentvertex==0) break;

        //******局部贪心策略:将当前未满足的节点currentvertex邻域中点度最大的一个非支配点变为支配点
		for(i=1;i<=gvertexnum;i++)
			if(GAdjMatrix[currentvertex][i] && DSList[i]==0)  //局部贪心地选择一个邻点(非支配点)--点度大优先(自然次序即表示度大的优先,这是因为图点序已按点度降序排了)
			{		bestvertex=i;
			        break;
			}
		if(i==gvertexnum+1)   //表示新添加任何邻点作为支配点都不会增加任何未满足的点的满足程度. 可以证明这种情况在不含孤立点的图中不会发生
		{
			printf("No solution to the instance!\n");
			return 0;
		}

		DSList[bestvertex]=1;          //加入新的支配点
		dssize++;                      //支配点个数更新
		for(j=1;j<=gvertexnum;j++)     //更新相关点的未满足程度Satisfied[.]
			if(GAdjMatrix[bestvertex][j])
			{    Satisfied[j]++;
  				 if(Satisfied[j]==0)   //如果点j被支配次数现在能达到其点度一半,则归于已满足的点类
					 uncoverednum--;
			}
	}


    //***用于调试阶段: 检查所得正影响的支配集是否正确
	if(!CheckPIDS(dssize))
	{	printf("There is some wrong with NewLocalGreedy3_PIDS!\n");
		getch();
	}//***用于调试阶段

	printf("\nNewLocalGreedy3_PIDS所求正影响支配集大小为(优化前): %5d\n",dssize);

    //优化阶段: 由所得正影响支配集导出一个极小正影响支配集
    int newdssize;
	newdssize=RefinePIDS(dssize);

	//***用于调试阶段: 检查所得正影响的支配集是否正确
	if(!CheckPIDS(newdssize))
	{	printf("There is some wrong with NewLocaGreedy3_PIDS!\n");
		getch();
	}//***用于调试阶段

	printf("\nNewLocalGreedy3_PIDS所求正影响支配集大小为(优化后): %5d\n",newdssize);
	return newdssize;//返回所求支配集大小
}







//***输入和输出函数***

void ReadGraph(char *txtgraphfilename)
{//To read the data file of graph G=<V,E>(V={1,2,...,n}) to obtain the adjacent matrix.
 //The format of the data file is as follows: the node number,the edge number, every edge <x,y> followed
      FILE  *in;

      //根据给定的文件名字打开文件。参数r表示对文件进行只读操作
	  if ((in=fopen(txtgraphfilename, "r"))==NULL)
      {  	fprintf(stderr, "Cannot open the data file.\n");
		return;
      }

      fscanf(in,"%d",&gvertexnum);                         //读取图的实际点数
      fscanf(in,"%d",&gedgenum);                           //读取图的边数

      int i,j,k;
      for(i=1;i<=gvertexnum;i++)                          //initialize the adjacent matrix
		for(j=1;j<=gvertexnum;j++)
			GAdjMatrix[i][j]=0;             /**/

      k=0;                                                 //k 用来记住实际的边数。注：文件中某边可能重复出现多次
      while(!feof(in))                                     //读出每一条边(每条边文件中是用一对点来表示的)
      {
	    fscanf(in,"%d",&i);
    	fscanf(in,"%d",&j);

	    if(GAdjMatrix[i][j])                              //如果边有了就不要重复。注：文件中边可能重复出现多次
             continue;
	    GAdjMatrix[i][j]=1;
	    GAdjMatrix[j][i]=GAdjMatrix[i][j];                //For an undirected graph, the adjacent matrix is symmetrical
	    k++;
      }
      gedgenum=k;                                         //得到实际的边数
      fclose(in);                                         //打开的文件最后必须关闭

     //初始化各点度,并处理孤立点(图中含有孤立点则问题无解). 检查是否含有孤立点. 如果含有孤立点,则将孤立点随机连边即可保证不含孤立点
  	 for(i=1;i<=gvertexnum;i++)
	 {  DegreeList[i]=0;
		for(j=1;j<=gvertexnum;j++)
			if(GAdjMatrix[i][j])
				DegreeList[i]++;
		if(DegreeList[i]==0)  //表示节点i是孤立点,此时将节点i随机连接到另一节点k
		{
           k=(rand()%gvertexnum)+1;
		   if(k==i && k<gvertexnum)
			   k++;
		   else if(k==i && k>1)
			   k--;
           GAdjMatrix[i][k]=1;
		   GAdjMatrix[k][i]=1;
		   DegreeList[i]++;
		   DegreeList[k]++;
		   gedgenum++;
		}
	 }

	 //对节点度按照由大到小排序, 然后对图按照度由大到小重新编号: 新1号节点度最大, 新最后第gvertexnum号节点度最小
	 //采用选择排序方法
	 short int temp;
	 bool tempf;
     int a;
	 for(i=1;i<gvertexnum;i++)
	 {
		 k=i;
		 for(j=i+1;j<=gvertexnum;j++)
			 if(DegreeList[j]>DegreeList[k])
			     k=j;
	     if(k!=i)
		 {//交换节点i和节点k的编号: 相应节点度要交换,同时要交换邻接矩阵的第i行和第k行,以及第i列和第k列
            temp=DegreeList[k];        //交换两节点的度
			DegreeList[k]=DegreeList[i];
			DegreeList[i]=temp;

			for(a=1;a<=gvertexnum;a++) //交换邻接矩阵的两行
			{	tempf=GAdjMatrix[k][a];
				GAdjMatrix[k][a]=GAdjMatrix[i][a];
				GAdjMatrix[i][a]=tempf;
			}

            for(a=1;a<=gvertexnum;a++) //交换邻接矩阵的两列
			{	tempf=GAdjMatrix[a][k];
				GAdjMatrix[a][k]=GAdjMatrix[a][i];
				GAdjMatrix[a][i]=tempf;
			}

		 }
	 }

	 //求最大度maxdegree, 最小度mindegree, 所有度为1的节点个数onedegreenum, 偶数度节点个数evendegreenum(奇数度节点个数为gvertexnum-evendegreenum)
	 maxdegree=DegreeList[1];
	 mindegree=DegreeList[gvertexnum];
	 onedegreenum=0;
	 evendegreenum=0;
	 for(i=gvertexnum;i>=1;i--)
		 if(DegreeList[i]==1)
			 onedegreenum++;
		 else if(DegreeList[i]%2==0)
			 evendegreenum++;

	 printf("n=%d,m=%d,maxd=%d,averd=%6.2lf,dens=%10.5lf\n",gvertexnum,gedgenum,maxdegree,2.0*gedgenum/gvertexnum,2.0*gedgenum/(gvertexnum*(gvertexnum-1.0)));
	 printf("mindegree=%d,onedegreenum=%d,evendegreenum=%d,odddegreenum=%d\n",mindegree,onedegreenum,evendegreenum,gvertexnum-evendegreenum);

}





//集成环境下执行部分(begin)****************************************************************
  int main(void)
  {

    char  graphfilename[50];
	char  resultfilename[50];
    //键盘上输入在ReadGraph()中要打开的图数据文件名(txt文件,带扩展名.文件格式:第一行图的点数n,第二行图的边数m,余下每一行是一条边(i,j)的一个节点i  j)
	printf("Please input the filename of the data file:\n");
	scanf("%s",graphfilename);
    //键盘上输入在下面中要创建的保存计算结果的文件名(txt文件,带扩展名)
	printf("Please input the filename of the result file:\n");
	scanf("%s",resultfilename);

//集成环境下执行部分(end)******************************************************************





	//------------------------------输入数据--------------------------------------

  	ReadGraph(graphfilename);  //读出图的数据
    //CreateGraph(graphfilename);  //随机产生一个图,并存储到文件中


	//-----------------------------执行算法--------------------------------------
	FILE *out;
	clock_t start, end;
    int dssize;

	//根据给定的结果文件名字打开文件。参数w表示对文件进行写操作
    if ((out=fopen(resultfilename, "w"))==NULL)
	{   	fprintf(stderr, "Cannot open the data file.\n");
			return 0;
	}
	fprintf(out,"******The graph file is %s******\n",graphfilename);
	fprintf(out,"vertexnum=%d,edgenum=%d,maxdegree=%d,averdegree=%6.2lf,density=%10.5lf\n",gvertexnum,gedgenum,maxdegree,2.0*gedgenum/gvertexnum,2.0*gedgenum/(gvertexnum*(gvertexnum-1.0)));
	fprintf(out, "mindegree=%d,onedegreenum=%d,evendegreenum=%d,odddegreenum=%d\n\n",mindegree,onedegreenum,evendegreenum,gvertexnum-evendegreenum);

    //******贪心近似算法1(文献中方法)*********
	printf("\nGreedy1_PIDS is running...");
    start=clock();
   	dssize=Greedy1_PIDS();      //执行贪心近似算法,返回所求支配集大小
	end=clock();
	if(dssize)
		fprintf(out,"\nGreedy1_PIDS obtains a PIDS of size %d.", dssize);
	else
		fprintf(out,"\nNo solution to the instance.");
	fprintf(out,"\nGreedy1_PIDS has running time %.3lf.\n", 1.0*(end-start)/CLK_TCK);
	//getch();

	/*
	//******贪心算法2(文献中方法)*********
	printf("\nGreedy2_PIDS is running...");
    start=clock();
   	dssize=Greedy2_PIDS();      //执行贪心近似算法,返回所求支配集大小
	end=clock();
	if(dssize)
		fprintf(out,"\nGreedy2_PIDS obtains a PIDS of size %d.", dssize);
	else
		fprintf(out,"\nNo solution to the instance.");
	fprintf(out,"\nGreedy2_PIDS has running time %.3lf.\n", 1.0*(end-start)/CLK_TCK);
	//getch();


	//*****新设计的贪心算法1*********
	printf("\nNewGreedy1_PIDS is running...");
    start=clock();
	dssize=NewGreedy1_PIDS();      //执行新设计的贪心算法1,返回所求支配集大小
	end=clock();
	if(dssize)
		fprintf(out,"\nNewGreedy1_PIDS obtains a PIDS of size %d.", dssize);
	else
		fprintf(out,"\nNo solution to the instance.");
    fprintf(out,"\nNewGreedy1_PIDS has running time %.3lf\n", 1.0*(end-start)/CLK_TCK);
    //getch();

	//*****新设计的贪心算法2*********
	printf("\nNewGreedy2_PIDS is running...");
    start=clock();
	dssize=NewGreedy2_PIDS();      //执行新设计的贪心算法2,返回所求支配集大小
	end=clock();
	if(dssize)
		fprintf(out,"\nNewGreedy2_PIDS obtains a PIDS of size %d.", dssize);
	else
		fprintf(out,"\nNo solution to the instance.");
    fprintf(out,"\nNewGreedy2_PIDS has running time %.3lf\n", 1.0*(end-start)/CLK_TCK);
    //getch();

	//*****新设计的贪心算法2*********
	printf("\nNewGreedy3_PIDS is running...");
    start=clock();
	dssize=NewGreedy3_PIDS();      //执行新设计的贪心算法2,返回所求支配集大小
	end=clock();
	if(dssize)
		fprintf(out,"\nNewGreedy3_PIDS obtains a PIDS of size %d.", dssize);
	else
		fprintf(out,"\nNo solution to the instance.");
    fprintf(out,"\nNewGreedy3_PIDS has running time %.3lf\n", 1.0*(end-start)/CLK_TCK);
    //getch();

	//*****新设计的局部贪心算法1*********
	printf("\nLocalGreedy1_PIDS is running...");
    start=clock();
	dssize=LocalGreedy1_PIDS();      //执行新设计的局部贪心算法1,返回所求支配集大小
	end=clock();
	if(dssize)
		fprintf(out,"\nLocalGreedy1_PIDS obtains a PIDS of size %d.", dssize);
	else
		fprintf(out,"\nNo solution to the instance.");
    fprintf(out,"\nLocalGreedy1_PIDS has running time %.3lf\n", 1.0*(end-start)/CLK_TCK);
    //getch();

	//*****新设计的局部贪心算法2*********
	printf("\nLocalGreedy2_PIDS is running...");
    start=clock();
	dssize=LocalGreedy2_PIDS();      //执行新设计的局部贪心算法2,返回所求支配集大小
	end=clock();
	if(dssize)
		fprintf(out,"\nLocalGreedy2_PIDS obtains a PIDS of size %d.", dssize);
	else
		fprintf(out,"\nNo solution to the instance.");
    fprintf(out,"\nLocalGreedy2_PIDS has running time %.3lf\n", 1.0*(end-start)/CLK_TCK);
    //getch();

	//*****新设计的局部贪心算法3*********
	printf("\nLocalGreedy3_PIDS is running...");
    start=clock();
	dssize=LocalGreedy3_PIDS();      //执行新设计的局部贪心算法3,返回所求支配集大小
	end=clock();
	if(dssize)
		fprintf(out,"\nLocalGreedy3_PIDS obtains a PIDS of size %d.", dssize);
	else
		fprintf(out,"\nNo solution to the instance.");
    fprintf(out,"\nLocalGreedy3_PIDS has running time %.3lf\n", 1.0*(end-start)/CLK_TCK);
    //getch();

	//*****新设计的新的局部贪心算法1*********
	printf("\nNewLocalGreedy1_PIDS is running...");
    start=clock();
	dssize=NewLocalGreedy1_PIDS();  //执行新设计的新的局部贪心算法1,返回所求支配集大小
	end=clock();
	if(dssize)
		fprintf(out,"\nNewLocalGreedy1_PIDS obtains a PIDS of size %d.", dssize);
	else
		fprintf(out,"\nNo solution to the instance.");
    fprintf(out,"\nNewLocalGreedy1_PIDS has running time %.3lf\n", 1.0*(end-start)/CLK_TCK);
    //getch();

	//*****新设计的新的局部贪心算法2*********
	printf("\nNewLocalGreedy2_PIDS is running...");
    start=clock();
	dssize=NewLocalGreedy2_PIDS();  //执行新设计的新的局部贪心算法2,返回所求支配集大小
	end=clock();
	if(dssize)
		fprintf(out,"\nNewLocalGreedy2_PIDS obtains a PIDS of size %d.", dssize);
	else
		fprintf(out,"\nNo solution to the instance.");
    fprintf(out,"\nNewLocalGreedy2_PIDS has running time %.3lf\n", 1.0*(end-start)/CLK_TCK);
    //getch();

	//*****新设计的新的局部贪心算法3*********
	printf("\nNewLocalGreedy3_PIDS is running...");
    start=clock();
	dssize=NewLocalGreedy3_PIDS();  //执行新设计的新的局部贪心算法3,返回所求支配集大小
	end=clock();
	if(dssize)
		fprintf(out,"\nNewLocalGreedy3_PIDS obtains a PIDS of size %d.", dssize);
	else
		fprintf(out,"\nNo solution to the instance.");
    fprintf(out,"\nNewLocalGreedy3_PIDS has running time %.3lf\n", 1.0*(end-start)/CLK_TCK);
    //getch();


    */

	fclose(out);           //打开的文件最后必须关闭
    printf("程序运行完毕!\n");
	//getch();
}































/*随机产生图的函数没有使用。啊哈哈哈*/
void CreateGraph(char *txtgraphfilename)
{//产生无向图,图中n个点1,2,...,n, 图的平均度averd=2m/n, 边的稠密度r=2m/(n(n-1)),因此有r=averd/(n-1).
  //图中任意点i和j之间存在边的概率为r(0<r<1)--r可以看作是边的稠密度.
   printf("Please input the number of nodes of the graph:\n");
   scanf("%d",&gvertexnum);       //确定图的点数

   double r;
   int averd;
   //printf("Please input the density r of edges of the graph (0<r<1):\n");
   //scanf("%lf",&r);             //确定图中边的稠密度
   printf("Please input the average degree of the graph (1<=averd<=n-1):\n");
	   scanf("%d", &averd);       //输入图的平均点度
   r=(1.0*averd)/(gvertexnum-1.0);           //确定图中边的稠密度

   int i,j;
   for(i=1;i<=gvertexnum;i++)     //初始化邻接矩阵
	   for(j=1;j<=gvertexnum;j++)
  	      GAdjMatrix[i][j]=0;

   time_t t;

   //随机产生图. 如果含有孤立点则将孤立点随机连边即可以保证不含孤立点

   srand((unsigned) time(&t));
   gedgenum=0;

   for(i=1;i<gvertexnum;i++)     //随机产生图的邻接矩阵
	   for(j=i+1;j<=gvertexnum;j++)
	   {
		   GAdjMatrix[i][j]=rand()%1000<(r*1000)?1:0;
	       if(GAdjMatrix[i][j])
		   {
			   GAdjMatrix[j][i]=1;
		       gedgenum++;
		   }
	   }

//初始化各点度,并处理孤立点(图中含有孤立点则问题无解). 检查是否含有孤立点. 如果含有孤立点,则将孤立点随机连边即可保证不含孤立点
	 int k;
  	 for(i=1;i<=gvertexnum;i++)
	 {  DegreeList[i]=0;
		for(j=1;j<=gvertexnum;j++)
			if(GAdjMatrix[i][j])
				DegreeList[i]++;
		if(DegreeList[i]==0)  //表示节点i是孤立点,此时将节点i随机连接到另一节点k
		{
           k=(rand()%gvertexnum)+1;
		   if(k==i && k<gvertexnum)
			   k++;
		   else if(k==i && k>1)
			   k--;
           GAdjMatrix[i][k]=1;
		   GAdjMatrix[k][i]=1;
		   DegreeList[i]++;
		   DegreeList[k]++;
		   gedgenum++;
		}
	 }

	 //对节点度按照由大到小排序, 然后对图按照度由大到小重新编号: 新1号节点度最大, 新最后第gvertexnum号节点度最小
	 //采用选择排序方法
	 short int temp;
	 bool tempf;
     int a;
	 for(i=1;i<gvertexnum;i++)
	 {
		 k=i;
		 for(j=i+1;j<=gvertexnum;j++)
			 if(DegreeList[j]>DegreeList[k])
			     k=j;
	     if(k!=i)
		 {//交换节点i和节点k的编号: 相应节点度要交换,同时要交换邻接矩阵的第i行和第k行,以及第i列和第k列
            temp=DegreeList[k];        //交换两节点的度
			DegreeList[k]=DegreeList[i];
			DegreeList[i]=temp;

			for(a=1;a<=gvertexnum;a++) //交换邻接矩阵的两行
			{	tempf=GAdjMatrix[k][a];
				GAdjMatrix[k][a]=GAdjMatrix[i][a];
				GAdjMatrix[i][a]=tempf;
			}

            for(a=1;a<=gvertexnum;a++) //交换邻接矩阵的两列
			{	tempf=GAdjMatrix[a][k];
				GAdjMatrix[a][k]=GAdjMatrix[a][i];
				GAdjMatrix[a][i]=tempf;
			}

		 }
	 }

	 //求最大度maxdegree, 最小度mindegree, 所有度为1的节点个数onedegreenum, 偶数度节点个数evendegreenum(奇数度节点个数为gvertexnum-evendegreenum)
	 maxdegree=DegreeList[1];
	 mindegree=DegreeList[gvertexnum];
	 onedegreenum=0;
	 evendegreenum=0;
	 for(i=gvertexnum;i>=1;i--)
		 if(DegreeList[i]==1)
			 onedegreenum++;
		 else if(DegreeList[i]%2==0)
			 evendegreenum++;

	 printf("n=%d,m=%d,maxd=%d,averd=%6.2lf,dens=%10.5lf\n",gvertexnum,gedgenum,maxdegree,2.0*gedgenum/gvertexnum,2.0*gedgenum/(gvertexnum*(gvertexnum-1.0)));
	 printf("mindegree=%d,onedegreenum=%d,evendegreenum=%d,odddegreenum=%d\n",mindegree,onedegreenum,evendegreenum,gvertexnum-evendegreenum);

	bool savegraph;
	savegraph=1;                      //savegraph=0 表示不用存储图到文件中
	if(savegraph)
	{	//存储图到指定的文件中
		FILE  *out;
		//根据函数传入的文件名字，创建并打开文件。参数w表示对文件进行写操作
		if ((out=fopen(txtgraphfilename, "w"))==NULL)
		{   	fprintf(stderr, "Cannot open the data file.\n");
			return;
		}
		fprintf(out,"%d\n",gvertexnum);
		fprintf(out,"%d\n",gedgenum);

		for(i=1;i<gvertexnum;i++)
			for(j=i+1;j<=gvertexnum;j++)
				if(GAdjMatrix[i][j])    //存储每条边
					fprintf(out,"%5d  %5d\n",i,j);
		fclose(out);                   //打开的文件最后必须关闭
	}
}


