$(DEST)/DMCCalculation.o: $(DEST)/ensemblewalkermodule.mod
$(DEST)/DMCCalculation.o: $(DEST)/inputreadermodule.mod
$(DEST)/DMCCalculation.o: $(DEST)/loggermodule.mod $(DEST)/numberkinds.mod
$(DEST)/DMCCalculation.o: $(DEST)/outputwritermodule.mod
$(DEST)/DMCCalculation.o: $(DEST)/potentialmodule.mod

$(DEST)/EnsembleWalker.o: $(DEST)/gaussiandistributormodule.mod
$(DEST)/EnsembleWalker.o: $(DEST)/histogrambuildermodule.mod
$(DEST)/EnsembleWalker.o: $(DEST)/loggermodule.mod $(DEST)/numberkinds.mod
$(DEST)/EnsembleWalker.o: $(DEST)/outputwritermodule.mod
$(DEST)/EnsembleWalker.o: $(DEST)/particleensemblemodule.mod
$(DEST)/EnsembleWalker.o: $(DEST)/potentialmodule.mod

$(DEST)/GaussianDistributor.o: $(DEST)/loggermodule.mod $(DEST)/numberkinds.mod
$(DEST)/GaussianDistributor.o: $(DEST)/outputwritermodule.mod

$(DEST)/GaussianDistributorTest.o: $(DEST)/gaussiandistributormodule.mod
$(DEST)/GaussianDistributorTest.o: $(DEST)/numberkinds.mod $(DEST)/testing.mod

$(DEST)/HistogramBuilder.o: $(DEST)/loggermodule.mod $(DEST)/numberkinds.mod
$(DEST)/HistogramBuilder.o: $(DEST)/outputwritermodule.mod

$(DEST)/HistogramBuilderTest.o: $(DEST)/histogrambuildermodule.mod
$(DEST)/HistogramBuilderTest.o: $(DEST)/numberkinds.mod $(DEST)/testing.mod

$(DEST)/InputReader.o: $(DEST)/loggermodule.mod $(DEST)/numberkinds.mod

$(DEST)/InputReaderTests.o: $(DEST)/inputreadermodule.mod
$(DEST)/InputReaderTests.o: $(DEST)/numberkinds.mod $(DEST)/testing.mod

$(DEST)/Logger.o: $(DEST)/numberkinds.mod $(DEST)/outputwritermodule.mod


$(DEST)/OutputWriter.o: $(DEST)/numberkinds.mod

$(DEST)/ParticleEnsemble.o: $(DEST)/loggermodule.mod $(DEST)/numberkinds.mod
$(DEST)/ParticleEnsemble.o: $(DEST)/outputwritermodule.mod

$(DEST)/ParticleEnsembleTests.o: $(DEST)/numberkinds.mod
$(DEST)/ParticleEnsembleTests.o: $(DEST)/particleensemblemodule.mod
$(DEST)/ParticleEnsembleTests.o: $(DEST)/testing.mod

$(DEST)/Potential.o: $(DEST)/loggermodule.mod $(DEST)/numberkinds.mod
$(DEST)/Potential.o: $(DEST)/outputwritermodule.mod

$(DEST)/PotentialTests.o: $(DEST)/numberkinds.mod $(DEST)/potentialmodule.mod
$(DEST)/PotentialTests.o: $(DEST)/testing.mod

$(DEST)/TestRunner.o: $(DEST)/gaussiandistributortests.mod
$(DEST)/TestRunner.o: $(DEST)/histogrambuildertests.mod
$(DEST)/TestRunner.o: $(DEST)/inputreadertests.mod $(DEST)/numberkinds.mod
$(DEST)/TestRunner.o: $(DEST)/particleensembletests.mod
$(DEST)/TestRunner.o: $(DEST)/potentialtests.mod

$(DEST)/Testing.o: $(DEST)/numberkinds.mod


$(DEST)/myapp.o: $(DEST)/dmccalculationmodule.mod

$(DEST)/dmccalculationmodule.mod: $(DEST)/DMCCalculation.o
$(DEST)/ensemblewalkermodule.mod: $(DEST)/EnsembleWalker.o
$(DEST)/gaussiandistributormodule.mod: $(DEST)/GaussianDistributor.o
$(DEST)/gaussiandistributortests.mod: $(DEST)/GaussianDistributorTest.o
$(DEST)/histogrambuildermodule.mod: $(DEST)/HistogramBuilder.o
$(DEST)/histogrambuildertests.mod: $(DEST)/HistogramBuilderTest.o
$(DEST)/inputreadermodule.mod: $(DEST)/InputReader.o
$(DEST)/inputreadertests.mod: $(DEST)/InputReaderTests.o
$(DEST)/loggermodule.mod: $(DEST)/Logger.o
$(DEST)/numberkinds.mod: $(DEST)/NumberKinds.o
$(DEST)/outputwritermodule.mod: $(DEST)/OutputWriter.o
$(DEST)/particleensemblemodule.mod: $(DEST)/ParticleEnsemble.o
$(DEST)/particleensembletests.mod: $(DEST)/ParticleEnsembleTests.o
$(DEST)/potentialmodule.mod: $(DEST)/Potential.o
$(DEST)/potentialtests.mod: $(DEST)/PotentialTests.o
$(DEST)/testing.mod: $(DEST)/Testing.o
$(DEST)/testrunnermodule.mod: $(DEST)/TestRunner.o
