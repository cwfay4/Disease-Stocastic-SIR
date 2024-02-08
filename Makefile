#
# Makefile for cre_prcltn
#
#
#all: mkgrph LoPR_ma Disease
all: mkgrph Disease Disease_m
 
COBJS=tri_prcltn.o graph.o cluster.o HK.o grphfns.o crprcfns.o clstrfns.o random.o
OBJS=mkgrph.o graph.o cluster.o HK.o random.o
LOBJS= LoPR_ma.o LoPR_alg_o.o random.o graph.o defaults.o LoPR_basic_utils.o
LoOBJS= LoPR.o LoPR_alg_o.o random.o graph.o defaults.o LoPR_basic_utils.o
DOBJS= Disease.o diseasePOP.o random.o graph.o Disease_basic_utils.o stats.o
DMOBJS= Disease_m.o diseasePOP.o random.o graph.o Disease_basic_utils.o stats.o

#CC=g++296
#CC=g++-4.9
#CC=g++ -std=c++0x
CC=g++
#CFLAGS = -g -o3 -Wall
#CFLAGS = -o3 -Wall
CFLAGS = -g -o3 
CXXFLAGS = -g
LIBS= -lm

tri_prcltn:  $(COBJS) $(COBJBIB)
	$(CC) $(CFLAGS) -o $@ $(COBJS) $(COBJBIB) $(LIBS)
mkgrph:  $(OBJS) $(OBJBIB)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(OBJBIB) $(LIBS)
LoPR_ma: 	 $(LOBJS) $(LOBJBIB)
	$(CC) $(CFLAGS) -o $@ $(LOBJS) $(LOBJBIB) $(LIBS) 
LoPR: 	 $(LoOBJS) $(LoOBJBIB)
	$(CC) $(CFLAGS) -o $@ $(LoOBJS) $(LoOBJBIB) $(LIBS) 
Disease: 	 $(DOBJS) $(DOBJBIB)
	$(CC) $(CFLAGS) -o $@ $(DOBJS) $(DOBJBIB) $(LIBS)
Disease_m: 	 $(DMOBJS) $(DMOBJBIB)
	$(CC) $(CFLAGS) -o $@ $(DMOBJS) $(DMOBJBIB) $(LIBS)	

clean:
	rm -f *.o tri_prcltn mkgrph LoPR_o_ma Disease Disease_m
