package util;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class FileUtils {
	
	public static List<String> fileLinesAsList(String file) throws IOException {
		List<String> rtrn = new ArrayList<String>();
		BufferedReader b = new BufferedReader(new FileReader(file));
		while(b.ready()) {
			rtrn.add(b.readLine());
		}
		b.close();
		return rtrn;
	}
	
}
