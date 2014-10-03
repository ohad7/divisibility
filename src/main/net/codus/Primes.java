package net.codus;
import java.io.File;
import java.math.BigInteger;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;


public class Primes {
	
	static private final List<BigInteger> primes = new ArrayList<>();
	
	static
	{
		try
		{
			System.out.println("Reading primes...");
			for (int i=1; i<=1; i++ ) {
				String filename = "src/main/resources/primes" + i +".txt";
				File file = new File(filename);
				System.out.print(" " + filename+".. ");
				List<String> lines = Files.readAllLines(file.toPath(), Charset.defaultCharset());
				for (String line : lines){
					for (String hey : line.split(" ")){
						if (hey.matches("\\d+"))
							primes.add(new BigInteger(hey.trim()));
					}
				}
			}
			System.out.println();
		}
		catch (Exception e) {
			throw new RuntimeException(e);
		}
	}
	
	public static BigInteger getPrime(int index){
		return primes.get(index);
	}
	

}
