package net.codus;
import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;


public class Benchmark {
	
	private static final int MAX_PRIME = 30000;
	
	public static void main(String[] args) throws IOException {
		
		
		Random random = new Random(0);
		Primes.getPrime(0);	// initialize the class
		
//		System.out.println("now press any key");
//		System.in.read();
//		Divisibility divisibility = new ModDivisibility(); // 1116ms
//		Divisibility divisibility = new BaseXDivisiblity();	// 18120ms
		Divisibility divisibility = new EfficientDivisibility(); // 35882ms // 25280ms (getLastDigit) //8957ms (divide)
		
		long startTime = System.currentTimeMillis();
		for (int i=0; i<100; i++){
			BigInteger number = BigInteger.ONE;
			List<BigInteger> primes = new ArrayList<>();
			for (int j=0; j<50; j++){
				BigInteger prime = Primes.getPrime(3+random.nextInt(MAX_PRIME));
				number = number.multiply(prime);
				primes.add(prime);
			}
			
			int primeIndex = 3;
			BigInteger prime = Primes.getPrime(primeIndex);
			System.out.println(i + ": number :" + number +" primes :" + primes);
			while (true){
//				if (i == 25 && prime.intValue() == 254491){
//					System.out.println(prime);
						
				boolean divisible = divisibility.isDivisibile(number, prime);
				//System.out.println(prime +" ->" + divisible);
				
//				if (divisible && !primes.contains(prime))
//					throw new RuntimeException("Should not be divisible :" + prime);
//				
//				if (!divisible && primes.contains(prime))
//					throw new RuntimeException("Should be divisible, but found not to be :" + prime);
//				}
				
				if (primeIndex++ > MAX_PRIME)
					break;
				prime = Primes.getPrime(primeIndex);
			}
			
			System.out.println("\n********************\n");
		}
		
		long endTime = System.currentTimeMillis();
		System.out.println("Time :" + (endTime-startTime)+"ms");
	}

}
