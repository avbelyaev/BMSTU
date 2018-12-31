import java.util.HashMap;

/**
 * Created by anthony on 09.10.16.
 */
public class VertexHashMap<K, V> extends HashMap<K, V> {

    @Override
    public V get(Object key) {
        V retval = null;
        K vertex = (K)key;

        super.forEach((k, v) -> {
            if (k.hashCode() == vertex.hashCode()){}
                //retval = super.get(vertex);
        });

        return retval;
    }
}
