package test;

import java.io.FileWriter;
import java.io.IOException;

import net.sf.samtools.util.CloseableIterator;
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;

public class TestUTRs {

	public static void main(String[] args) throws IOException {
		
		AnnotationCollection<Gene> genes = BEDFileIO.loadFromFile(args[0], args[1]);
		CloseableIterator<Gene> iter = genes.sortedIterator();
		
		FileWriter w3 = new FileWriter(args[0] + "_3UTR.bed");
		FileWriter w5 = new FileWriter(args[0] + "_5UTR.bed");
		
		while(iter.hasNext()) {
			Gene gene = iter.next();
			try {
				Annotation utr5 = gene.get5UTR();
				Annotation utr3 = gene.get3UTR();
				String line3 = utr3.toBED() + "\n";
				String line5 = utr5.toBED() + "\n";
				w3.write(line3);
				w5.write(line5);
				if(line3.contains("null")) {
					System.out.println(gene.getName() + "\tUTR3");
				}
				if(line5.contains("null")) {
					System.out.println(gene.getName() + "\tUTR5");
				}
			} catch(NullPointerException e) {
				continue;
			}
		}
		
		iter.close();
		w3.close();
		w5.close();
		
	}
	
}
