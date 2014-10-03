package net.codus;
import java.math.BigInteger;


public class ModDivisibility implements Divisibility{

	@Override
	public boolean isDivisibile(BigInteger number, BigInteger prime) {
		return number.mod(prime).signum() == 0;
	}

}
