package file;

import general.CommandLineParser;
import general.StringParser;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;


/**
 * Write the text of a Java method to write out the contents of a provided file
 * So the file can be written from within code instead of being passed as input to program
 * @author prussell
 *
 */
public class HardCodedFileWriter {
	
	public static void writeMethodBodyToFile(String originalFile, String outFile) throws IOException {
		FileReader r = new FileReader(originalFile);
		BufferedReader b = new BufferedReader(r);
		FileWriter w = new FileWriter(outFile);
		StringParser s = new StringParser();
		while(b.ready()) {
			String line = b.readLine();
			s.parse(line);
			String toWrite = "w.write(\"";
			for(int i = 0; i < s.getFieldCount() - 1; i++) {
				toWrite += s.asString(i) + "\\t";
			}
			toWrite += s.asString(s.getFieldCount() - 1);
			toWrite += "\\n\");";
			w.write(toWrite + "\n");
		}
		r.close();
		b.close();
		w.close();
	}
	
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("--orig", "Original file", true);
		p.addStringArg("--out", "File to write method body to", true);
		p.parse(args);
		String orig = p.getStringArg("--orig");
		String out = p.getStringArg("--out");
		
		writeMethodBodyToFile(orig, out);
		
	}
	

}
